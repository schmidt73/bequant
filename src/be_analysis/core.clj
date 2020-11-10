(ns be-analysis.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cheshire.core :as chesire]))

;;;; Data loading code 

(defn get-sensor
  [s]
  (let [lscaffold "AAAAAGTGGCACCGAGTCGGTGCTTTTTTT"
        rscaffold "GAATTC"]
    (if-let [lindex (string/last-index-of s lscaffold)]
      (if-let [rindex (string/last-index-of s rscaffold)]
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

(defn load-guide-csv-row
  [row]
  (let [edit-pos (Integer/parseInt (nth row 18))
        target (subs (nth row 13) 135 175)]
    {:sg-rna (nth row 12)
     :guide-id (nth row 7)
     :target target
     :edit-pos edit-pos}))

(defn lazily-load-guides
  [rdr]
  (->> (rest (csv/read-csv rdr))
       (map load-guide-csv-row)))

;;;; Data analysis code

(defn get-outcomes
  [guide sensor]
  (->> (map #(if (not= %1 %2) {:from %1 :to %2}) (:target guide) sensor)
       (map-indexed #(if %2 (assoc %2 :position %1)))
       (filter identity)))

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
        (assoc-in outcomes-map [guide :total] (inc total-count))
        outcomes)))
   {}
   guide-sensor-seq))

;;;; Pretty printing code

(defn pretty-print-outcomes-edit-pos
  [outcomes writer]
  (.append writer "Guide ID,Sequence,C-to-A,C-to-T,C-to-G,Total\n")
  (doseq [guide (keys outcomes)]
    (let [outcome (get outcomes guide)
          guide-id (:guide-id guide)
          edit-pos (+ (:edit-pos guide) 9)
          total (:total outcome 0)
          sequence (:target guide)
          c-to-a (get outcome {:from \C :to \A :position edit-pos} 0)
          c-to-t (get outcome {:from \C :to \T :position edit-pos} 0)
          c-to-g (get outcome {:from \C :to \G :position edit-pos} 0)
          line (string/join "," [guide-id sequence c-to-a c-to-t c-to-g total])]
      (.append writer (str line "\n")))))

(defn jsonify-outcomes
  [outcomes]
  (map (fn [[guide outcomes]]
         [guide
          {:total (:total outcomes)
           :outcomes (->> (map (fn [[k v]] (if (map? k) (assoc k :count v))) outcomes)
                          (filter identity))}])
       outcomes))

(defn pretty-print-outcomes-json
  [outcomes writer]
  (-> (jsonify-outcomes outcomes)
      (chesire/generate-stream writer)))

;;;; Code that actually does data analysis

(defn report-count-filter [msg f coll]
  (let [result (filter f coll)]
    (println (str msg (count result)))
    result))

;; (if-let [lindex (string/last-index-of s lscaffold)]
;;   (if-let [rindex (string/last-index-of s rscaffold)]
;;     (if (< lindex rindex) 
;;       (subs s (+ lindex (count lscaffold)) rindex)))))

(defn get-sg-rna
  [s]
  (let [lscaffold "CACC"
        rscaffold "GTTTAAG"]
    (if-let [lindex (string/index-of s lscaffold)]
      (if-let [rindex (string/index-of s rscaffold)]
        (if (< lindex rindex)
          (subs s (+ lindex (count lscaffold)) rindex))))))

(defn count-statistics-sensor-nn-matching
  [guides fastq-file]
  "Counts various statistics for a read file."
  (let [sg-rna-to-guide (into {} (map #(vector (:sg-rna %) %) guides))
        lsh (load-lsh guides)]
    (with-open [reader (io/reader fastq-file)]
      (->> (line-seq reader)
           (map-indexed (fn [idx item] (if (= 0 (mod (dec idx) 4)) item)))
           (report-count-filter "# reads: " identity)
           (map get-sensor)
           (report-count-filter "# reads with sensors: " some?)
           (report-count-filter "# reads with length 40 sensors: "
                                #(= 40 (count %)))
           (pmap #(find-approx-nearest-neighbor lsh % 3))
           (report-count-filter "# reads with length 40 sensors and matched guides: "
                                some?)
           (doall)))))

(defn count-statistics-sgrna-matching
  [guides fastq-file]
  "Counts various statistics for a read file."
  (let [sg-rna-to-guide (into {} (map #(vector (:sg-rna %) %) guides))
        lsh (load-lsh guides)]
    (with-open [reader (io/reader fastq-file)]
      (->> (line-seq reader)
           (map-indexed (fn [idx item] (if (= 0 (mod (dec idx) 4)) item)))
           (report-count-filter "# reads: " identity)
           (map #(vector (get sg-rna-to-guide (get-sg-rna %))
                         (get-sensor %)))
           (report-count-filter "# reads with sensors: "
                                (fn [[_ sensor]] (some? sensor)))
           (report-count-filter "# reads with length 40 sensors: "
                                (fn [[_ sensor]] (= 40 (count sensor))))
           (report-count-filter "# reads with length 40 sensors and matched guides: "
                                (fn [[guides _]] (some? guides)))
           (doall)))))

(defn analyze-fastq-file-with-progress-exact-sgrna
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
           (filter #(= 40 (count (second %))))
           (map #(do
                   (when (= 0 (mod @counter 10000))
                     (println (str "Progress: " @counter)))
                   (swap! counter inc)
                   %))
           (count-outcomes)))))

(defn analyze-fastq-file-with-progress
  [lsh guides fastq-file]
  (let [counter (atom 0)]
    (with-open [reader (io/reader fastq-file)]
      (->> (lazily-load-sensors reader)
           (map #(do
                   (when (= 0 (mod @counter 10000))
                     (println (str "Progress: " @counter)))
                   (swap! counter inc)
                   %))
           (analyze-sensors lsh guides)))))

(def parent-dir
  "/home/schmidt73/Desktop/base-editing/")

(def hbes-merged-fastq-files
  (->> (io/file (str parent-dir "/mdamb231-fastq/merged/"))
       (file-seq)
       (filter #(.isFile %))
       (filter #(re-find #"HBES" (.getName %)))))

(def mbes-merged-fastq-files
  (->> (io/file (str parent-dir "/mdamb231-fastq/merged/"))
       (file-seq)
       (filter #(.isFile %))
       (filter #(re-find #"MBES" (.getName %)))))

(def hbes-guides-csv-file
  (str parent-dir "/mdamb231-fastq/HBESv4_revised_whitelist.csv"))

(def mbes-guides-csv-file
  (str parent-dir "/mdamb231-fastq/MBESv4_revised_whitelist.csv"))

(defn load-guides-file
  [guides-csv-file]
  (with-open [rdr (io/reader guides-csv-file)]
    (->> (lazily-load-guides rdr)
         (vec))))

(defn load-lsh
  [guides]
  (create-lsh-structure 25 10 40 (map :target guides))) 

(-main)

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
