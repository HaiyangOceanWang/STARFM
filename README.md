# STARFM

(1) Filtering out spectrally different pixels and cloud pixels in the moving window.

Two ways to filter out spectrally different pixels:
•	unsupervised classification
•	thresholds in surface reflectance directly
Cloud-free pixels

(2) Filtering out pixels those couldn’t provide better spectral and spatial information than central pixel.
•	All poor-quality data are excluded from candidates according to the QA layer in the Landsat and MODIS surface reflectance products.
•	refine

(3) Both fine- and coarse-resolution data are then used in determining sample weights (Wijk).

(4) Compute the prediction value using 
