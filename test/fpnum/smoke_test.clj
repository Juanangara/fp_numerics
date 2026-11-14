(ns fpnum.smoke-test
  "Prueba de humo: confirma que el entorno de pruebas funciona y que todos los
  namespaces del proyecto se cargan sin errores."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.test.check :as tc]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]))

(def namespaces
  '[fpnum.linear
    fpnum.ad.protocol
    fpnum.ad.forward
    fpnum.ad.ops
    fpnum.ad.reverse
    fpnum.ad.syntax
    fpnum.svd.core
    fpnum.svd.image
    fpnum.optim
    fpnum.demo.ad
    fpnum.demo.conditioning
    fpnum.demo.descent
    fpnum.demo.reverse
    fpnum.demo.svd])

(deftest namespaces-se-cargan
  (doseq [ns-sym namespaces]
    (testing (str "carga " ns-sym)
      (is (do (require ns-sym) true)))))

(deftest demos-tienen-main
  (doseq [ns-sym '[fpnum.demo.ad fpnum.demo.conditioning fpnum.demo.descent fpnum.demo.reverse fpnum.demo.svd]]
    (require ns-sym)
    (is (some? (ns-resolve ns-sym '-main)) (str ns-sym " debe definir -main"))))

(deftest test-check-funciona
  (testing "propiedad trivial: invertir dos veces un vector devuelve el original"
    (let [propiedad (prop/for-all [v (gen/vector gen/small-integer)]
                      (= v (vec (reverse (reverse v)))))
          resultado (tc/quick-check 100 propiedad)]
      (is (:pass? resultado)))))
