(ns be-analysis.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cheshire.core :as chesire]))

(defn hamming-distance
  [s1 s2]
  (count (filter not (map = s1 s2))))

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
    {:sg-rna (nth row 11) ; 12 for HBES 11 for MBES
     :guide-id (nth row 7)
     :target target
     :pam (nth row 16)
     :edit-pos edit-pos}))

(defn load-guide-csv-hbes-row
  [row]
  (let [edit-pos (Integer/parseInt (nth row 18))
        target (subs (nth row 13) 135 175)]
    {:sg-rna (nth row 12) ; 12 for HBES 11 for MBES
     :guide-id (nth row 7)
     :target target
     :pam (nth row 15)
     :edit-pos edit-pos}))

(defn lazily-load-guides
  [rdr row-reader]
  (->> (rest (csv/read-csv rdr))
       (map row-reader)))

;;;; Data analysis code

(defn get-outcomes
  [guide sensor]
  (if (= 40 (count sensor))
    (->> (map #(if (not= %1 %2) {:from %1 :to %2}) (:target guide) sensor)
         (map-indexed #(if %2 (assoc %2 :position %1)))
         (filter identity)
         (#(conj %2 %1) :edits {:type :edit :distance (hamming-distance guide sensor)}))
    [:indels {:type :indel :size (- (count sensor) 40)}]))

(defn count-outcomes
  [guide-sensor-seq]
  (reduce
   (fn [outcomes-map [guide sensor]]
     (let [outcomes (get-outcomes guide sensor)
           total-count (get-in outcomes-map [guide :total] 0)]
       (reduce
        (fn [counts outcome]
          (let [old-count (get-in outcomes-map [guide outcome] 0)]
            (assoc-in counts [guide outcome] (inc old-count))))
        outcomes-map
        outcomes)))
   {}
   guide-sensor-seq))

;;;; Pretty printing code

(defn pretty-print-row
  [guide outcome-rep1 outcome-rep2 writer]
  (let [guide-id (:guide-id guide)
        edit-pos (+ (:edit-pos guide) 9)
        sequence (:target guide)
        total-rep1  (+ (get outcome-rep1 :edits 0) (get outcome-rep1 :indels 0))
        percent-rep1 #(float (if (= total-rep1 0) 0 (* 100 (/ % total-rep1))))
        c-to-a-rep1 (get outcome-rep1 {:from \C :to \A :position edit-pos} 0)
        c-to-t-rep1 (get outcome-rep1 {:from \C :to \T :position edit-pos} 0)
        c-to-g-rep1 (get outcome-rep1 {:from \C :to \G :position edit-pos} 0)
        indels-rep1 (get outcome-rep1 :indels 0)
        total-rep2  (+ (get outcome-rep2 :edits 0) (get outcome-rep2 :indels 0))
        percent-rep2 #(float (if (= total-rep2 0) 0 (* 100 (/ % total-rep2))))
        c-to-a-rep2 (get outcome-rep2 {:from \C :to \A :position edit-pos} 0)
        c-to-t-rep2 (get outcome-rep2 {:from \C :to \T :position edit-pos} 0)
        c-to-g-rep2 (get outcome-rep2 {:from \C :to \G :position edit-pos} 0)
        indels-rep2 (get outcome-rep2 :indels 0)
        line (string/join "," [guide-id total-rep1 c-to-t-rep1 c-to-a-rep1 c-to-g-rep1
                               indels-rep1 total-rep2 c-to-t-rep2 c-to-a-rep2 c-to-g-rep2
                               indels-rep2 (percent-rep1 c-to-t-rep1) (percent-rep1 c-to-a-rep1)
                               (percent-rep1 c-to-g-rep1) (percent-rep1 indels-rep1)
                               (percent-rep2 c-to-t-rep2) (percent-rep2 c-to-a-rep2)
                               (percent-rep2 c-to-g-rep2) (percent-rep2 indels-rep2)])]
    (.append writer (str line "\n"))))

(defn pretty-print-outcomes-edit-pos-row
  [guide outcome-rep1 outcome-rep2 edit-pos writer]
  (let [guide-id (:guide-id guide)
        pam (:pam guide)
        nuc-context (subs (:target guide) (max (- edit-pos 2) 0) (+ edit-pos 3))
        total-rep1  (+ (get outcome-rep1 :edits 0) (get outcome-rep1 :indels 0))
        percent-rep1 #(float (if (= total-rep1 0) 0 (* 100 (/ % total-rep1))))
        c-to-a-rep1 (get outcome-rep1 {:from \C :to \A :position edit-pos} 0)
        c-to-t-rep1 (get outcome-rep1 {:from \C :to \T :position edit-pos} 0)
        c-to-g-rep1 (get outcome-rep1 {:from \C :to \G :position edit-pos} 0)
        total-rep2  (+ (get outcome-rep2 :edits 0) (get outcome-rep2 :indels 0))
        percent-rep2 #(float (if (= total-rep2 0) 0 (* 100 (/ % total-rep2))))
        c-to-a-rep2 (get outcome-rep2 {:from \C :to \A :position edit-pos} 0)
        c-to-t-rep2 (get outcome-rep2 {:from \C :to \T :position edit-pos} 0)
        c-to-g-rep2 (get outcome-rep2 {:from \C :to \G :position edit-pos} 0)
        line (string/join "," [guide-id edit-pos pam nuc-context
                               total-rep1 c-to-t-rep1 c-to-a-rep1 c-to-g-rep1
                               total-rep2 c-to-t-rep2 c-to-a-rep2 c-to-g-rep2
                               (percent-rep1 c-to-t-rep1) (percent-rep1 c-to-a-rep1)
                               (percent-rep1 c-to-g-rep1) (percent-rep2 c-to-t-rep2)
                               (percent-rep2 c-to-a-rep2) (percent-rep2 c-to-g-rep2)])]
    (.append writer (str line "\n"))))

(defn pretty-print-edit-pos-rows
  [guide outcome-rep1 outcome-rep2 writer]
  (doseq [edit-pos (->> (map-indexed #(vector %1 %2) (:target guide))
                        (filter #(and (<= (first %) 30)
                                      (= (second %) \C)))
                        (map first))]
    (pretty-print-outcomes-edit-pos-row
     guide outcome-rep1 outcome-rep2 edit-pos writer)))

(defn pretty-print-outcomes-csv
  [outcomes-rep1 outcomes-rep2 writer]
  (let [header (str "guide_ID,total_REP1,tCTN_REP1,tCAN_REP1,tCGN_REP1,"
                    "indel_REP1,total_REP2,tCTN_REP2,tCAN_REP2,tCGN_REP2,"
                    "indel_REP2,percent_tCTN_REP1,percent_tCAN_REP1,"
                    "percent_tCGN_REP1,percent_indel_REP1,"
                    "percent_tCTN_REP2,percent_tCAN_REP2,percent_tCGN_REP2,"
                    "percent_indel_REP2,percent_tCTN,percent_tCAN,percent_tCGN,"
                    "percent_indel")]
    (.append writer (str header "\n"))
    (doseq [guide (keys outcomes-rep1)]
      (let [outcome-rep1 (get outcomes-rep1 guide)
            outcome-rep2 (get outcomes-rep2 guide)]
        (pretty-print-row guide outcome-rep1 outcome-rep2 writer)))))

(defn pretty-print-outcomes-edit-pos-csv
  [outcomes-rep1 outcomes-rep2 writer]
  (let [header (str "guide_ID,cytosine_position,PAM,surrounding nucleotide context (NNCNN),"
                    "total_REP1,tCTN_REP1,tCAN_REP1,tCGN_REP1,total_REP2,tCTN_REP2,tCAN_REP2,"
                    "tCGN_REP2,percent_tCTN_REP1,percent_tCAN_REP1,percent_tCGN_REP1,"
                    "percent_tCTN_REP2,percent_tCAN_REP2,percent_tCGN_REP2,percent_tCTN,"
                    "percent_tCAN,percent_tCGN")]
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
                     (println (str "Progress: " @counter)))
                   (swap! counter inc)
                   %))
           (count-outcomes)))))

(doseq [base-editor (keys hbes-merged-fastq-files)]
  (let [reps (get hbes-merged-fastq-files base-editor)
        rep1 (first (filter #(re-find #"REP1" (.getName %)) reps))
        rep2 (first (filter #(re-find #"REP2" (.getName %)) reps))
        guides (load-hbes-guides-file hbes-guides-csv-file)
        outcomes-rep1 (analyze-fastq-file-with-progress guides rep1)
        outcomes-rep2 (analyze-fastq-file-with-progress guides rep2)
        output-file (str (.getName rep1) ".csv")]
    (with-open [writer (io/writer output-file)]
      (pretty-print-outcomes-edit-pos-csv outcomes-rep1 outcomes-rep2 writer))))

(def parent-dir
  "/home/schmidt73/Desktop/base-editing/")

(def hbes-merged-fastq-files
  (->> (io/file (str parent-dir "/mdamb231-fastq/merged/"))
       (file-seq)
       (filter #(.isFile %))
       (filter #(re-find #"HBES" (.getName %)))
       (group-by #(second (re-find #"HBES_(.*)_REP" (.getName %))))))
       
(def mbes-merged-fastq-files
  (->> (io/file (str parent-dir "/mdamb231-fastq/merged/"))
       (file-seq)
       (filter #(.isFile %))
       (filter #(re-find #"MBES" (.getName %)))
       (group-by #(second (re-find #"MBES_(.*)_REP" (.getName %))))))

(def hbes-guides-csv-file
  (str parent-dir "/mdamb231-fastq/HBESv4_revised_whitelist.csv"))

(def mbes-guides-csv-file
  (str parent-dir "/mdamb231-fastq/MBESv4_revised_whitelist.csv"))

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

(defn load-lsh
  [guides]
  (create-lsh-structure 25 10 40 (map :target guides))) 

(defn -main [& args]
  (doseq [mbes-fastq mbes-merged-fastq-files]
    (with-open [reader (io/reader mbes-fastq)]
      (do
        (println (str "Analyzing: " mbes-fastq))
        (let [guides (load-guides-file mbes-guides-csv-file)
              outcomes (analyze-fastq-file-with-progress-exact-sgrna guides mbes-fastq)
              output-file (str (.getName mbes-fastq) ".json")]
          (println (str "Writing outcomes to file: " output-file))
          (with-open [writer (io/writer output-file)]
            (pretty-print-outcomes-json outcomes writer))))))
  (doseq [hbes-fastq hbes-merged-fastq-files]
    (with-open [reader (io/reader hbes-fastq)]
      (do
        (println (str "Analyzing: " hbes-fastq "\n"))
        (let [guides (load-guides-file hbes-guides-csv-file)
              outcomes (analyze-fastq-file-with-progress-exact-sgrna guides hbes-fastq)
              output-file (str (.getName hbes-fastq) ".json")]
          (println (str "Writing outcomes to file: " output-file))
          (with-open [writer (io/writer output-file)]
            (pretty-print-outcomes-json outcomes writer)))))))
