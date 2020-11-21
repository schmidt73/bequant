(ns be-analysis.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [taoensso.timbre :as timbre]
            [cheshire.core :as chesire]
            [clojure.tools.cli :refer [parse-opts]]))

;;;; Data loading code 

(defn get-sensor
  [s]
  (let [lscaffold "AAAAAGTGGCACCGAGTCGGTGCTTTTTTT"
        rscaffold "GAATTC"]
    (if-let [lindex (string/last-index-of s lscaffold)]
      (if-let [rindex (string/last-index-of s rscaffold)]
        (if (< lindex rindex) 
          (subs s (+ lindex (count lscaffold)) rindex))))))

(defn get-sg-rna
  [s]
  (let [lscaffold "CACC"
        rscaffold "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT"]
    (if-let [lindex (string/index-of s lscaffold)]
      (if-let [rindex (string/index-of s rscaffold)]
        (if (< lindex rindex)
          (subs s (+ lindex (count lscaffold)) rindex))))))

(defn lazily-load-sensors
  [rdr]
  (->> (line-seq rdr)
       (map-indexed (fn [idx item] (if (= 0 (mod (dec idx) 4)) item)))
       (filter identity)
       (map get-sensor)
       (filter identity)
       (filter #(= 40 (count %)))))

(defn load-guide-csv-mbes-row
  [row]
  (let [edit-pos (Integer/parseInt (nth row 18))
        target (subs (nth row 13) 135 175)]
    {:sg-rna (nth row 11)
     :guide-id (nth row 7)
     :target target
     :pam (nth row 31)
     :edit-pos edit-pos}))

(defn load-guide-csv-hbes-row
  [row]
  (let [edit-pos (Integer/parseInt (nth row 18))
        target (subs (nth row 13) 135 175)]
    {:sg-rna (nth row 12)
     :guide-id (nth row 7)
     :target target
     :pam (nth row 31)
     :edit-pos edit-pos}))

(defn lazily-load-guides
  [rdr row-reader]
  (->> (rest (csv/read-csv rdr))
       (map row-reader)))

(defn load-mbes-guides-file
  [guides-csv-file]
  (with-open [rdr (io/reader guides-csv-file)]
    (->> (lazily-load-guides rdr load-guide-csv-mbes-row)
         (vec))))

(defn load-hbes-guides-file
  [guides-csv-file]
  (with-open [rdr (io/reader guides-csv-file)]
    (->> (lazily-load-guides rdr load-guide-csv-hbes-row)
         (vec))))

;;;; Data analysis code

(defn hamming-distance
  [s1 s2]
  (count (filter not (map = s1 s2))))

(defn get-outcomes
  [guide sensor]
  (let [hd (hamming-distance (:target guide) sensor)]
    (cond
        (not= 40 (count sensor))
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

(defn pretty-print-row
  [guide outcome-rep1 outcome-rep2 writer]
  (let [guide-id (:guide-id guide)
        edit-pos (+ (:edit-pos guide) 9)
        sequence (:target guide)
        total-rep1  (+ (get outcome-rep1 :edits 0) (get outcome-rep1 :indels 0)
                       (get outcome-rep1 :non-edits 0))
        percent-rep1 #(float (if (= total-rep1 0) 0 (* 100 (/ % total-rep1))))
        c-to-a-rep1 (total-edit-count outcome-rep1 \C \A edit-pos)
        c-to-t-rep1 (total-edit-count outcome-rep1 \C \T edit-pos)
        c-to-g-rep1 (total-edit-count outcome-rep1 \C \G edit-pos)
        c-to-t-nc-rep1 (get outcome-rep1 {:from \C :to \T :position edit-pos :distance 1} 0)
        indels-rep1 (get outcome-rep1 :indels 0)
        total-rep2  (+ (get outcome-rep2 :edits 0) (get outcome-rep2 :indels 0)
                       (get outcome-rep2 :non-edits 0))
        percent-rep2 #(float (if (= total-rep2 0) 0 (* 100 (/ % total-rep2))))
        c-to-a-rep2 (total-edit-count outcome-rep2 \C \A edit-pos)
        c-to-t-rep2 (total-edit-count outcome-rep2 \C \T edit-pos)
        c-to-g-rep2 (total-edit-count outcome-rep2 \C \G edit-pos)
        c-to-t-nc-rep2 (get outcome-rep2 {:from \C :to \T :position edit-pos :distance 1} 0)
        indels-rep2 (get outcome-rep2 :indels 0)
        line (string/join "," [guide-id
                               total-rep1 c-to-t-rep1 c-to-a-rep1 c-to-g-rep1
                               c-to-t-nc-rep1 indels-rep1
                               total-rep2 c-to-t-rep2 c-to-a-rep2 c-to-g-rep2
                               c-to-t-nc-rep2 indels-rep2
                               (percent-rep1 c-to-t-rep1) (percent-rep1 c-to-a-rep1)
                               (percent-rep1 c-to-g-rep1) (percent-rep1 c-to-t-nc-rep1)
                               (percent-rep1 indels-rep1)
                               (percent-rep2 c-to-t-rep2) (percent-rep2 c-to-a-rep2)
                               (percent-rep2 c-to-g-rep2) (percent-rep2 c-to-t-nc-rep2)
                               (percent-rep2 indels-rep2)])]
    (.append writer (str line "\n"))))

(defn pretty-print-outcomes-edit-pos-row
  [guide outcome-rep1 outcome-rep2 edit-pos writer]
  (let [guide-id (:guide-id guide)
        pam (:pam guide)
        nuc-context (subs (:target guide) (max (- edit-pos 2) 0) (+ edit-pos 3))
        total-rep1  (+ (get outcome-rep1 :edits 0) (get outcome-rep1 :indels 0))
        percent-rep1 #(float (if (= total-rep1 0) 0 (* 100 (/ % total-rep1))))
        c-to-a-rep1 (total-edit-count outcome-rep1 \C \A edit-pos)
        c-to-t-rep1 (total-edit-count outcome-rep1 \C \T edit-pos)
        c-to-g-rep1 (total-edit-count outcome-rep1 \C \G edit-pos)
        c-to-t-nc-rep1 (get outcome-rep1 {:from \C :to \T :position edit-pos :distance 1} 0)
        total-rep2  (+ (get outcome-rep2 :edits 0) (get outcome-rep2 :indels 0))
        percent-rep2 #(float (if (= total-rep2 0) 0 (* 100 (/ % total-rep2))))
        c-to-a-rep2 (total-edit-count outcome-rep2 \C \A edit-pos)
        c-to-t-nc-rep2 (get outcome-rep2 {:from \C :to \T :position edit-pos :distance 1} 0)
        c-to-t-rep2 (total-edit-count outcome-rep2 \C \T edit-pos)
        c-to-g-rep2 (total-edit-count outcome-rep2 \C \G edit-pos)
        line (string/join "," [guide-id edit-pos pam nuc-context
                               total-rep1 c-to-t-rep1 c-to-a-rep1 c-to-g-rep1 c-to-t-nc-rep1
                               total-rep2 c-to-t-rep2 c-to-a-rep2 c-to-g-rep2 c-to-t-nc-rep2
                               (percent-rep1 c-to-t-rep1) (percent-rep1 c-to-a-rep1)
                               (percent-rep1 c-to-g-rep1) (percent-rep1 c-to-t-nc-rep1)
                               (percent-rep2 c-to-t-rep2) (percent-rep2 c-to-a-rep2)
                               (percent-rep2 c-to-g-rep2) (percent-rep2 c-to-t-nc-rep2)])]
    (.append writer (str line "\n"))))

(defn pretty-print-edit-pos-rows
  [guide outcome-rep1 outcome-rep2 writer]
  (doseq [edit-pos (->> (map-indexed #(vector %1 %2) (:target guide))
                        (filter #(and (<= (first %) 30)
                                      (>= (first %) 5)
                                      (= (second %) \C)))
                        (map first))]
    (pretty-print-outcomes-edit-pos-row
     guide outcome-rep1 outcome-rep2 edit-pos writer)))

(defn pretty-print-outcomes-csv
  [outcomes-rep1 outcomes-rep2 writer]
  (let [header (str "guide_ID,"
                    "total_REP1,tCTN_REP1,tCAN_REP1,tCGN_REP1,tCT_REP1,indel_REP1,"
                    "total_REP2,tCTN_REP2,tCAN_REP2,tCGN_REP2,tCT_REP2,indel_REP2,"
                    "percent_tCTN_REP1,percent_tCAN_REP1,"
                    "percent_tCGN_REP1,percent_tCT_REP1,percent_indel_REP1,"
                    "percent_tCTN_REP2,percent_tCAN_REP2,"
                    "percent_tCGN_REP2,percent_tCT_REP2,percent_indel_REP2")]
    (.append writer (str header "\n"))
    (doseq [guide (keys outcomes-rep1)]
      (let [outcome-rep1 (get outcomes-rep1 guide)
            outcome-rep2 (get outcomes-rep2 guide)]
        (pretty-print-row guide outcome-rep1 outcome-rep2 writer)))))

(defn pretty-print-outcomes-edit-pos-csv
  [outcomes-rep1 outcomes-rep2 writer]
  (let [header (str "guide_ID,cytosine_position,PAM,surrounding nucleotide context (NNCNN),"
                    "total_REP1,tCTN_REP1,tCAN_REP1,tCGN_REP1,tCT_REP1,"
                    "total_REP2,tCTN_REP2,tCAN_REP2,tCGN_REP2,tCT_REP2,"
                    "percent_tCTN_REP1,percent_tCAN_REP1,percent_tCGN_REP1,percent_tCT_REP1,"
                    "percent_tCTN_REP2,percent_tCAN_REP2,percent_tCGN_REP2,percent_tCT_REP2")]
    (.append writer (str header "\n"))
    (doseq [guide (keys outcomes-rep1)]
      (let [outcome-rep1 (get outcomes-rep1 guide)
            outcome-rep2 (get outcomes-rep2 guide)]
        (pretty-print-edit-pos-rows guide outcome-rep1 outcome-rep2 writer)))))

(defn jsonify-outcomes
  [outcomes]
  (map (fn [[guide outcomes]]
         [guide
          {:total (+ (get outcomes :edits 0) (get outcomes :indels 0))
           :indels (get outcomes :indels 0)
           :outcomes (->> (map (fn [[k v]] (if (map? k) (assoc k :count v))) outcomes)
                          (filter identity))}])
       outcomes))

(defn pretty-print-outcomes-json
  [outcomes writer]
  (-> (jsonify-outcomes outcomes)
      (chesire/generate-stream writer)))

;;;; Code that actually does data analysis

(defn analyze-fastq-file-with-progress
  [guides fastq-file]
  (let [counter (atom 0)
        sg-rna-to-guide (into {} (map #(vector (:sg-rna %) %) guides))]
    (with-open [reader (io/reader fastq-file)]
      (->> (line-seq reader)
           (map-indexed (fn [idx item] (if (= 0 (mod (dec idx) 4)) item)))
           (filter identity)
           (map #(vector (get sg-rna-to-guide (get-sg-rna %))
                         (get-sensor %)))
           (filter #(and (some? (first %)) (some? (second %))))
           (map #(do
                   (when (= 0 (mod @counter 10000))
                     (timbre/info (str "FASTq Records Processed: " @counter)))
                   (swap! counter inc)
                   %))
           (count-outcomes)))))

(def guides-map
  {:mbes (load-mbes-guides-file (io/resource "MBESv4_revised_whitelist.csv"))
   :hbes (load-hbes-guides-file (io/resource "HBESv4_revised_whitelist.csv"))})

(defn process-fastqs-edit-position
  [output-file rep1 rep2 screen-type format]
  (let [guides (screen-type guides-map)
        outcomes-rep1 (do (timbre/info "Processing FASTq file:" rep1)
                          (analyze-fastq-file-with-progress guides rep1))
        outcomes-rep2 (do (timbre/info "Processing FASTq file:" rep2)
                          (analyze-fastq-file-with-progress guides rep2))]
    (with-open [writer (io/writer output-file)]
      (case format
        :target (pretty-print-outcomes-csv outcomes-rep1 outcomes-rep2 writer)
        :all (pretty-print-outcomes-edit-pos-csv outcomes-rep1 outcomes-rep2 writer)))))
  
;;;; CLI 

(def cli-options
  [["-o" "--output FILE" "Output file (REQUIRED)."]
   ["-s" "--screen SCREEN" "Screen to pull whitelist for"
    :default :mbes
    :default-desc "MBES"
    :parse-fn #(keyword (string/lower-case %))
    :validate [#{:mbes :hbes} "must be one of ['MBES', 'HBES']"]]
   ["-f" "--format OUTPUT_FORMAT" "Output either only target cytosine or all cytosines."
    :default :target
    :default-desc "TARGET"
    :parse-fn #(keyword (string/lower-case %))
    :validate [#{:target :all} "must be one of ['TARGET', 'ALL']"]]
   ["-h" "--help"]])

(defn usage [options-summary]
  (->> ["Usage: java -jar analyze-fastqs rep1 rep2 [options]"
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

      (and (= 2 (count arguments)) (:output options))
      {:action [(:output options) (first arguments) (second arguments)
                (:screen options) (:format options)]}

      :else {:exit-message (usage summary)})))

(defn exit [status msg]
  (println msg)
  (System/exit status)) 

(defn -main [& args]
  (let [{:keys [action options exit-message ok?]} (validate-args args)]
    (if exit-message
      (exit (if ok? 0 1) exit-message)
      (apply process-fastqs-edit-position action))))
