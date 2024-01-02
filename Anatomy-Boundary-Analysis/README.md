# SizePelvis

## Overview

This MATLAB project provides advanced methods for comparing the size of pelvic structures using multiple methodologies. It includes boundary boxes, boundary spheres, and landmark distances to provide a comprehensive analysis of the size and proportions of pelvic structures in different datasets. It can also be transferred to other anatomical structures and models.

## Key Features:
1. **Multiple Size Methods**: This code includes several size comparison methods:
- Bounding Box: Two different approaches are used for the bounding box comparison — the hull method and the Singular Value Decomposition (SVD) method.
- Boundary Sphere: The code estimates the size of the pelvis using a boundary sphere encapsulating the structure.
- Landmark Distances: As a reference (or ground truth), the distances between specific anatomical landmarks are used. An alternate method, termed "ground truth 2," uses landmark distances to the acetabulum centre.
2. **Size and Scale Calculations**: Each method provides metrics to compute the size and scaling factor for each pelvis data set.
3. **Correlation Analysis**: The program evaluates which landmarks determine the size and scaling factor the most, based on Pearson correlation coefficients.
4. **Deviation Analysis**: An analysis is performed on the variance or standard deviation of size and scale factors to compare the methods.

## How to Use:
1. **Code Structure**: The main script is designed to control and access the different classes and functions. Key parts of the code are extracted into individual functions for further use and reusability.
2. **Setup** Ensure you have MATLAB installed. Clone this project to your local machine.
3. **Input Data**: Place your pelvis data files in the designated \GeometriesPelves directory. The main script expects data in a specific format (stl).
4. **Run**: Launch MATLAB and navigate to the project directory. Run the main script.
5. **Results**: The results, representing comparisons across all methodologies, will be displayed in MATLAB and saved in the \Figures directory.

## Use Cases:
- Ideal for biomedical engineering, particularly orthopaedic research.
- Useful for scaling and comparing anatomical models in medical studies.

## Benefits:
- Provides multiple scaling methods for improved analysis accuracy.
- Efficient and comprehensive visualisation and analysis of anatomical models.

## Compatibility:
The class is compatible with MATLAB and can be easily integrated into existing workflows for scaling anatomical models.

## Disclaimer:
This MATLAB class is provided on the MATLAB File Exchange for educational and research purposes. Users should ensure that the class meets their specific analysis requirements and may need to adapt it accordingly. The code is provided "as-is," and the author assumes no responsibility for its use or any consequences thereof.

## MATLAB Version Compatibility:
Optimised for Matlab 2023a, but should be compatible with most recent MATLAB versions.

## Keywords:
Pelvis, Landmarks, Size Comparison, Scaling, Bounding Box, Bounding Sphere, Biomedical Modeling, STL, SVD, convhull.

## References:
- [Scaling Methods of the Pelvis](https://doi.org/10.1515/cdbme-2022-1203) - Micheler, C. M. et al. (2022)
- [Minimal Bounding Box](https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box) - Johannes Korsawe (2024)
- [calc_OriBoundingBox(data)](https://www.mathworks.com/matlabcentral/fileexchange/64417-calc_oriboundingbox-data) - Svenja Reimer (2024)
- [Exact minimum bounding spheres and circles](https://github.com/AntonSemechko/Bounding-Spheres-And-Circles) - Anton Semechko (2024)
