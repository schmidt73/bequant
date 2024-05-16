(ns be-analysis.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [clojure.data.json :as json]
            [clojure.java.io :as io]
            [taoensso.timbre :as timbre]
            [cheshire.core :as chesire]
            [clojure.math.combinatorics :as combinatorics]
            [clojure.tools.cli :refer [parse-opts]]))

;;;; Data loading code 
(defn get-sensor
  [config s]
  (let [lscaffold (:sensor-left-scaffold config)
        rscaffold (:sensor-right-scaffold config)] 
    (if-let [lindex (string/last-index-of s lscaffold)]
      (if-let [rindex (string/last-index-of s rscaffold)]
        (if (< lindex rindex) 
          (subs s (+ lindex (count lscaffold)) rindex))))))

(defn get-sg-rna
  [config s]
  (let [lscaffold (:sgrna-left-scaffold config)
        rscaffold (:sgrna-right-scaffold config)]
    (if-let [lindex (string/index-of s lscaffold)]
      (if-let [rindex (string/index-of s rscaffold)]
        (if (< lindex rindex)
          (subs s (+ lindex (count lscaffold)) rindex))))))

(defn load-guide-csv-row
  [config base row]
  (let [edit-pos (Integer/parseInt (nth row 18))
        target (subs (nth row 13) (:target-start config) (:target-end config))]
    {:sg-rna (nth row 12)
     :guide-id (nth row 7)
     :target target
     :pam (nth row 31)
     :edit-pos edit-pos}))

(defn lazily-load-guides
  [rdr row-reader]
  (->> (rest (csv/read-csv rdr))
       (map row-reader)))

(defn load-whitelist
  [row-reader guides-csv-file]
  (with-open [rdr (io/reader guides-csv-file)]
    (->> (lazily-load-guides rdr row-reader)
         (vec))))

(defn load-config
  [config-file]
  (-> (slurp config-file)
      (json/read-str :key-fn keyword)))

(defn jsonify-outcomes
  [outcomes]
  (into {}
    (map (fn [[guide outcomes]]
           [(:guide-id guide)
            {:total (:total outcomes)
             :outcomes (->> (map (fn [[k v]] (if (map? k) (assoc k :count v))) outcomes)
                            (filter identity))}])
     outcomes)))

(defn pretty-print-outcomes-json
  [outcomes writer]
  (-> (jsonify-outcomes outcomes)
      (chesire/generate-stream writer)))

;;;; Data analysis code
(defn hamming-distance
  [s1 s2]
  (count (filter not (map = s1 s2))))

;; TODO: FIX, it incorrectly detects INDELS at moment
(defn get-outcomes
  [guide sensor]
  (let [hd (hamming-distance (:target guide) sensor)]
    (cond
        (not= (count (:target guide)) (count sensor))
        [:indels]

        (= 0 hd)
        [:non-edits]

        (not= 0 (hamming-distance (subs sensor 0 5) (subs (:target guide) 0 5)))
        [:guide-sensor-mismatches]

        :otherwise
        (->> (map #(if (not= %1 %2) {:from %1 :to %2 :distance hd}) (:target guide) sensor)
             (map-indexed #(if %2 (assoc %2 :position %1)))
             (filter identity)
             (concat [:edits])))))

(defn count-outcomes
  [guide-sensor-seq]
  (reduce
   (fn [outcomes-map [guide sensor]]
     (let [outcomes (get-outcomes guide sensor)]
       (reduce
        (fn [counts outcome]
          (let [old-count (get-in outcomes-map [guide outcome] 0)]
            (assoc-in counts [guide outcome] (inc old-count))))
        outcomes-map
        outcomes)))
   {}
   guide-sensor-seq))

;;;; Pretty printing code

(defn total-edit-count
  [outcome from to position]
  (reduce #(+ %1 (get outcome {:from from :to to :position position :distance %2} 0))
          0 (range 45)))

(defn pretty-print-outcome-row
  [edit-pos base outcome]
  (let [total  (+ (get outcome :edits 0)
                  (get outcome :indels 0)
                  (get outcome :non-edits 0))
        non-edits "N/A"
        base-to-a (if (= \A base) non-edits (total-edit-count outcome base \A edit-pos))
        base-to-c (if (= \C base) non-edits (total-edit-count outcome base \C edit-pos))
        base-to-g (if (= \G base) non-edits (total-edit-count outcome base \G edit-pos))
        base-to-t (if (= \T base) non-edits (total-edit-count outcome base \T edit-pos))
        base-to-a-nc (if (= \A base) non-edits (get outcome {:from base :to \A :position edit-pos :distance 1} 0))
        base-to-c-nc (if (= \C base) non-edits (get outcome {:from base :to \C :position edit-pos :distance 1} 0))
        base-to-g-nc (if (= \G base) non-edits (get outcome {:from base :to \G :position edit-pos :distance 1} 0)) 
        base-to-t-nc (if (= \T base) non-edits (get outcome {:from base :to \T :position edit-pos :distance 1} 0))
        indels (get outcome :indels 0)]
    (string/join "," [total base-to-a base-to-c base-to-g base-to-t
                      base-to-a-nc base-to-c-nc base-to-g-nc base-to-t-nc indels])))

(defn pretty-print-row
  [guide base outcomes writer]
  (let [guide-id (:guide-id guide)
        edit-pos (+ (:edit-pos guide) 9)
        sequence (:target guide)
        printed-outcomes (map #(pretty-print-outcome-row edit-pos base (get % guide)) outcomes)
        line (str guide-id "," (string/join "," printed-outcomes))]
    (.append writer (str line "\n"))))

(defn pretty-print-outcome-edit-pos-row
  [outcome base edit-pos]
  (let [total  (+ (get outcome :edits 0) (get outcome :indels 0))
        non-edits "N/A"
        base-to-a (if (= \A base) non-edits (total-edit-count outcome base \A edit-pos))
        base-to-c (if (= \C base) non-edits (total-edit-count outcome base \C edit-pos))
        base-to-g (if (= \G base) non-edits (total-edit-count outcome base \G edit-pos))
        base-to-t (if (= \T base) non-edits (total-edit-count outcome base \T edit-pos))
        base-to-a-nc (if (= \A base) non-edits (get outcome {:from base :to \A :position edit-pos :distance 1} 0))
        base-to-c-nc (if (= \C base) non-edits (get outcome {:from base :to \C :position edit-pos :distance 1} 0))
        base-to-g-nc (if (= \G base) non-edits (get outcome {:from base :to \G :position edit-pos :distance 1} 0))
        base-to-t-nc (if (= \T base) non-edits (get outcome {:from base :to \T :position edit-pos :distance 1} 0))
        indels (get outcome :indels 0)]
    (string/join "," [total base-to-a base-to-c base-to-g base-to-t
                      base-to-a-nc base-to-c-nc base-to-g-nc base-to-t-nc indels])))

(defn pretty-print-outcomes-edit-pos-row
  [guide base outcomes edit-pos writer]
  (let [guide-id (:guide-id guide)
        sequence (:sg-rna guide)
        pam (:pam guide)
        ex-pam (subs (:target guide) (+ 31 (count pam)) (+ 33 (count pam)))
        nuc-context (subs (:target guide) (max (- edit-pos 2) 0) (+ edit-pos 3))
        prefix (string/join "," [guide-id sequence pam ex-pam edit-pos nuc-context])
        suffix (string/join "," (map #(pretty-print-outcome-edit-pos-row (get % guide) base edit-pos) outcomes))]
    (.append writer (str prefix "," suffix "\n"))))

(defn pretty-print-edit-pos-rows
  [guide base outcomes writer]
  (doseq [edit-pos (->> (map-indexed #(vector %1 %2) (:target guide))
                        (filter #(and (<= (first %) 30)
                                      (>= (first %) 5)
                                      (= (second %) base)))
                        (map first))]
    (pretty-print-outcomes-edit-pos-row guide base outcomes edit-pos writer)))

(def pp-outcomes-header-format-string
    {\C (str "total_REP%d,tCAN_REP%d,tCCN_REP%d,tCGN_REP%d,tCTN_REP%d,"
             "tCA_REP%d,tCC_REP%d,tCG_REP%d,tCT_REP%d,tINDEL_REP%d"),
     \T (str "total_REP%d,tTAN_REP%d,tTCN_REP%d,tTGN_REP%d,tTTN_REP%d,"
             "tTA_REP%d,tTC_REP%d,tTG_REP%d,tTT_REP%d,tINDEL_REP%d"),
     \G (str "total_REP%d,tGAN_REP%d,tGCN_REP%d,tGGN_REP%d,tGTN_REP%d,"
             "tGA_REP%d,tGC_REP%d,tGG_REP%d,tGT_REP%d,tINDEL_REP%d")
     \A (str "total_REP%d,tAAN_REP%d,tACN_REP%d,tAGN_REP%d,tATN_REP%d,"
             "tAA_REP%d,tAC_REP%d,tAG_REP%d,tAT_REP%d,tINDEL_REP%d")})

(defn pretty-print-outcomes-csv
  [outcomes writer base]
  (let [outcome-headers
        (for [i (range (count outcomes))]
          (apply format
                 (get pp-outcomes-header-format-string base)
                 (take 1000 (repeat (+ 1 i)))))
        header (str "guide_ID," (string/join "," outcome-headers))]
    (.append writer (str header "\n"))
    (doseq [guide (keys (first outcomes))]
      (pretty-print-row guide base outcomes writer))))

(defn pretty-print-outcomes-edit-pos-csv
  [outcomes writer base]
  (let [outcome-headers
        (for [i (range (count outcomes))]
          (apply format
                 (get pp-outcomes-header-format-string base)
                 (take 1000 (repeat i))))
        header (str "guide_ID,sequence,PAM,exPAM,"
                    "edit_position,surrounding nucleotide context (NNCNN),"
                    (string/join "," outcome-headers))]
    (.append writer (str header "\n"))
    (doseq [guide (keys (first outcomes))]
      (pretty-print-edit-pos-rows guide base outcomes writer))))

(defn display [x] (prn x) (prn (count x)) x)

;; Given a sgRNA and a base, return all possible sgRNA 
;; sequences that can be generated by changing the base.
(defn expand-sg-rna
  [sg-rna base]
  (let [indices (filterv #(not (nil? %)) (map-indexed (fn [idx c] (if (= c base) idx)) sg-rna))
        nucleotides [\A \C \G \T]
        possibilities (map #(into {} (map vector indices %)) (apply combinatorics/cartesian-product (repeat (count indices) nucleotides)))]
      (when (< 7 (count indices))
        (timbre/warn (str "Warning: sgRNA suffix " sg-rna " has " (count indices) " " base " nucleotides. This will result in a large number of possible sgRNAs.")))
      (for [possibility possibilities]
        (string/join (map-indexed (fn [idx c] (if (contains? possibility idx) (get possibility idx) c)) sg-rna)))))

;rewrite this idiomatically using threading macros (into {} (apply concat )))

(defn expand-sg-rnas-in-whitelist
  [guides base]
  (let [counter (atom 0)
        total-guides (count guides)
        process-guide (fn [guide]
                        (when (= 0 (mod @counter 10))
                          (timbre/info (str "Processed " @counter "/" total-guides " guides in whitelist.")))
                        (swap! counter inc)
                        (map vector (expand-sg-rna (:sg-rna guide) base) (repeat guide)))
        expanded-sg-rnas (map process-guide guides)]
    (into {} (apply concat expanded-sg-rnas))))

(defn analyze-fastq-file-with-progress
  [config guides base fastq-file]
  (timbre/info "Expanding sgRNA sequences in whitelist...")
  (let [counter (atom 0)
        sg-rna-to-guide (expand-sg-rnas-in-whitelist guides base)] ; a dictionary that takes in the sgrna and returns the matching guide
    (timbre/info "Analyzing FASTq file:" fastq-file)
    (with-open [reader (io/reader fastq-file)]
      (->> (line-seq reader)
           (map-indexed (fn [idx item] (if (= 0 (mod (dec idx) 4)) item)))
           (filter identity)
           (map #(vector (get sg-rna-to-guide (get-sg-rna config %))
                         (get-sg-rna config %)
                         (get-sensor config %)))
           (filter #(and (some? (first %)) (some? (second %))))
           (map #(do
                   (when (= 0 (mod @counter 10000))
                     (timbre/info (str "FASTq Records Processed: " @counter)))
                   (swap! counter inc)
                   %))
           (count-outcomes)))))

(defn process-fastqs-edit-position
  [output-file replicates whitelist format config base]
  (let [config (load-config config)
        guides (load-whitelist #(load-guide-csv-row config base %) whitelist)
        outcomes (for [rep replicates]
                   (do (timbre/info "Processing FASTq file:" rep)
                       (analyze-fastq-file-with-progress config guides base rep)))]
    (with-open [writer (io/writer output-file)]
      (case format
        :target (pretty-print-outcomes-csv outcomes writer base)
        :json (pretty-print-outcomes-json outcomes writer)
        :all (pretty-print-outcomes-edit-pos-csv outcomes writer base)))))
  
;;;; CLI 
(def cli-options
  [["-o" "--output FILE" "Output file (REQUIRED)."]
   ["-w" "--whitelist WHITELIST" "Whitelist file for screen (REQUIRED)."]
   ["-c" "--config CONFIG" "Configuration file describing sensor and whitelist structure (REQUIRED)."]
   ["-b" "--base BASE" "Base to analyze (REQUIRED)."
    :default \C
    :default-desc "C"
    :parse-fn #(get % 0)
    :validate [#{\A \C \T \G} "must be one of ['A', 'T', 'C', 'G']"]]
   ["-f" "--format OUTPUT_FORMAT" "Output either only target cytosine or all cytosines."
    :default :target
    :default-desc "TARGET"
    :parse-fn #(keyword (string/lower-case %))
    :validate [#{:target :all :json} "must be one of ['TARGET', 'ALL', 'JSON']"]]
   ["-h" "--help"]])

(defn usage [options-summary]
  (->> ["Usage: java -jar analyze-fastqs rep1 ... repN [options]"
        ""
        "Options:"
        options-summary
        ""
        "Expects FASTq files consisting of the paired end reads"
        "for each replicate."]
       (string/join \newline)))

(defn error-msg [errors]
  (str "The following errors occurred while parsing your command:\n\n"
       (string/join \newline errors)))

(defn validate-args
  [args]
  (let [{:keys [options arguments summary errors]} (parse-opts args cli-options)]
    (cond
      (:help options) {:exit-message (usage summary) :ok? true}
      errors          {:exit-message (error-msg errors)}
      (and (<= 1 (count arguments)) (:output options) (:whitelist options) (:config options))
      {:action [(:output options) arguments 
                (:whitelist options) (:format options) (:config options) (:base options)]}

      :else {:exit-message (usage summary)})))

(defn exit [status msg]
  (println msg)
  (System/exit status)) 

(defn -main [& args]
  (let [{:keys [action options exit-message ok?]} (validate-args args)]
    (if exit-message
      (exit (if ok? 0 1) exit-message)
      (apply process-fastqs-edit-position action))))
