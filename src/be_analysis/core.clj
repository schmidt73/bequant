(ns be-analysis.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io])
  (:import java.security.MessageDigest))

;;; LSH code

(defn hamming-distance
  [s1 s2]
  (count (filter not (map = s1 s2))))

(defn hash-func
  "Returns a hash function that uses CRC-32 over the sampled bits to
  compute a hash."
  [bits-to-sample]
  (fn [msg]
    (let [crc (new java.util.zip.CRC32)]
      (doseq [n bits-to-sample]
        (.update crc (int (nth msg n))))
      (.getValue crc))))

(defn create-lsh-structure
  "Creates an LSH structure for NN search with l hash
  functions that sample k bits in d dimensions."
  [l k d strings]
  (let [create-hash-func (fn [] (hash-func (take k (repeatedly #(rand-int d)))))
        hash-funcs (vec (take l (repeatedly create-hash-func)))
        create-table (fn [hash-func] (into {} (map #(vector (hash-func %1) %2) strings strings)))
        hash-tables (mapv create-table hash-funcs)]
    {:hash-funcs hash-funcs
     :hash-tables hash-tables}))

(defn search-lsh
  [lsh idx string]
  (let [hash-func (nth (:hash-funcs lsh) idx)
        hash-table (nth (:hash-tables lsh) idx)]
    (if-let [match (get hash-table (hash-func string))]
      {:distance (hamming-distance match string)
       :match match})))

(defn find-approx-nearest-neighbor
  [lsh string r]
  (->> (map #(search-lsh lsh % string) (range (count (:hash-funcs lsh))))
       (filter identity)
       (some #(if (<= (:distance %) r) %))))

;;;; Data loading code 

(def lscaffold "AAAAAGTGGCACCGAGTCGGTGCTTTTTTT")
(def rscaffold "GAATTC")

(defn get-sensor
  [s]
  (if-let [lindex (string/last-index-of s lscaffold)]
    (if-let [rindex (string/last-index-of s rscaffold)]
      (if (< lindex rindex) 
        (subs s (+ lindex (count lscaffold)) rindex)))))

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
    {:guide-id (nth row 7)
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

(defn analyze-sensors
  "Analyzes a set of sensors, returning a map from the guides
  to the number of edits in each position."
  [lsh guides sensors]
  (let [targets-to-guides (into {} (map #(vector (:target %) %)) guides)]
     (->> sensors
          (pmap #(vector % (find-approx-nearest-neighbor lsh % 3)))
          (filter second)
          (pmap #(vector (get targets-to-guides (:match (second %))) (first %)))
          (count-outcomes))))

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

;;;; Code that actually does data analysis

;; (let [guides (load-guides-file hbes-guides-csv-file)
;;       lsh (load-lsh guides)]
;;   (pretty-print-outcomes-edit-pos
;;    (analyze-fastq-file-with-progress lsh guides (first hbes-merged-fastq-files))
;;    *out*))

(defn analyze-fastq-file-with-progress
  [lsh guides fastq-file]
  (let [counter (atom 0)]
    (with-open [reader (io/reader fastq-file)]
      (->> (lazily-load-sensors reader)
           (map #(do
                   (when (= 0 (mod @counter 10000))
                     (prn (str "Progress: " @counter)))
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

(defn -main [& args]
  (doseq [mbes-fastq mbes-merged-fastq-files]
    (with-open [reader (io/reader mbes-fastq)]
      (do
        (print (str "Analyzing: " mbes-fastq))
        (let [guides (load-guides-file mbes-guides-csv-file)
              lsh (load-lsh guides)
              outcomes (analyze-fastq-file-with-progress mbes-fastq)
              output-file (str (.getName mbes-fastq) "csv")]
          (print (str "Writing outcomes to file: " mbes-fastq))
          (with-open [writer (io/writer output-file)]
            (pretty-print-outcomes-edit-pos outcomes writer))))))
  (doseq [hbes-fastq hbes-merged-fastq-files]
    (with-open [reader (io/reader hbes-fastq)]
      (do
        (print (str "Analyzing: " hbes-fastq))
        (let [guides (load-guides-file hbes-guides-csv-file)
              lsh (load-lsh guides)
              outcomes (analyze-fastq-file-with-progress hbes-fastq)
              output-file (str (.getName hbes-fastq) "csv")]
          (print (str "Writing outcomes to file: " hbes-fastq))
          (with-open [writer (io/writer output-file)]
            (pretty-print-outcomes-edit-pos outcomes writer)))))))
            
      
