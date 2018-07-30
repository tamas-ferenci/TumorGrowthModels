# TumorGrowthMixedEffects
Simple mixed-effects growth models applied to tumor growth (modeling the effect of treatment and
measurement method).

These models - or their variants - were used in the following papers:
* Johanna Sápi, Tamás Ferenci, Dániel András Drexler, Levente Kovács. Tumor model identification
and statistical analysis. In: IEEE International Conference on Systems, Man, and Cybernetics 2015:
IEEE SMC 2015. Hong Kong, China, 2015.10.09-2015.10.12. Hong Kong: IEEE, 2015. pp. 2481-2486.
(ISBN:978-1-4799-8697-2).
* T Ferenci, J Sápi, L Kovács. Modelling xenograft tumor growth under antiangiogenic inhibitation
with mixed-effects models. In: 2016 IEEE International Conference on Systems, Man, and Cybernetics
Conference Proceedings: SMC 2016. Budapest, Hungary, 2016.10.09-2016.10.12. Budapest: IEEE, 2016.
pp. 3912-3917. (ISBN:978-1-5090-1897-0).
* Tamás Ferenci, Johanna Sápi, Levente Kovács. Modelling Tumor Growth Under Angiogenesis Inhibition
with Mixed-effects Models. ACTA POLYTECHNICA HUNGARICA 14:(1) pp. 221-234. (2017).

Code is available [here](https://github.com/tamas-ferenci/TumorGrowthMixedEffects/blob/master/TumorGrowthMixedEffects.R).
Columns of `RawData` are: `Date` (date of the measurement), `Type` (`C` denotes control, `E` denotes
treated), `Code` (unique identifier of the mouse), `Caliper1`, `Caliper2`, `Caliper3`, `MRI`
(measurements of the given mouse at the given date with the respective measurement methods).