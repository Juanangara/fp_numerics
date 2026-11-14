(ns fpnum.optim-grad-fn-test
  "Pruebas de la opción :grad-fn de fpnum.optim/gradient-descent con el gradiente
  de forward (por defecto) y el de reverse.

  Viven aparte de optim_test.clj porque dependen de fpnum.ad.reverse; así cada
  commit del historial puede ejecutar sus pruebas sin necesitar el módulo reverse."
  (:require [clojure.test :refer [deftest is testing]]
            [fpnum.ad.forward :as ad]
            [fpnum.ad.reverse :as rev]
            [fpnum.linear :as la]
            [fpnum.optim :as opt]))

(defn- approx-num= [a b tol]
  (<= (Math/abs (- (double a) (double b)))
      (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- approx-vec= [x y tol]
  (and (= (count x) (count y))
       (every? true? (map #(approx-num= %1 %2 tol) x y))))
(deftest grad-fn-por-defecto-es-forward
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 100 :w-true [2.0 -1.0 0.5] :sigma 0.1 :seed 42})
        f (opt/least-squares-objective x y)
        opts {:alpha (opt/safe-step x) :max-iter 5000 :tol 1e-8}
        sin-opcion (opt/gradient-descent f (la/zeros 3) opts)
        explicito (opt/gradient-descent f (la/zeros 3) (assoc opts :grad-fn ad/grad))]
    (is (= sin-opcion explicito) "omitir :grad-fn equivale a pasar fpnum.ad.forward/grad")))

(deftest descenso-con-reverse-da-lo-mismo-que-con-forward
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 100 :w-true [2.0 -1.0 0.5] :sigma 0.1 :seed 42})
        f (opt/least-squares-objective x y)
        opts {:alpha (opt/safe-step x) :max-iter 5000 :tol 1e-8}
        con-forward (opt/gradient-descent f (la/zeros 3) opts)
        con-reverse (opt/gradient-descent f (la/zeros 3) (assoc opts :grad-fn rev/grad))]
    (is (:converged? con-reverse))
    (testing "mismas iteraciones y mismo punto (los gradientes coinciden a ~1e-12)"
      (is (= (:iterations con-forward) (:iterations con-reverse)))
      (is (approx-vec= (:x con-forward) (:x con-reverse) 1e-9)))
    (testing "y la misma solución que las ecuaciones normales"
      (is (approx-vec= (opt/normal-equations x y) (:x con-reverse) 1e-6)))))

(deftest descenso-con-backtracking-y-reverse-sobre-rosenbrock
  (let [opts {:backtracking? true :max-iter 300 :tol 1e-6}
        con-forward (opt/gradient-descent opt/rosenbrock [-1.2 1.0] opts)
        con-reverse (opt/gradient-descent opt/rosenbrock [-1.2 1.0] (assoc opts :grad-fn rev/grad))]
    (is (= (:iterations con-forward) (:iterations con-reverse)))
    (is (approx-vec= (:x con-forward) (:x con-reverse) 1e-9))
    (is (every? true? (map (fn [a b] (approx-num= (:alpha a) (:alpha b) 1e-9))
                           (butlast (:trajectory con-forward))
                           (butlast (:trajectory con-reverse))))
        "el paso elegido por backtracking en cada iteración también coincide")))
