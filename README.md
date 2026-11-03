# HAT-Coast-GEE
Selected Google Earth Engine code for the HAT-Coast framework

## Associated manuscript

HAT-Coast: A Multi-source hierarchical adaptive thresholding framework for biennial global mainland shoreline mapping from 2000 to 2020

## Scope

This repository provides selected methodological functions used for:

- locally adaptive optical and SAR water extraction;
- hierarchical monthly gap filling;
- annual Water Inundation Frequency calculation; and
- persistent-water classification using the validation-derived threshold.

## Input conventions

For locally adaptive optical and SAR masks:

- 0 = non-water
- 1 = water

For the final persistent-water map:

- 0 = land or intermittently inundated area
- 1 = persistent water

## Citation

Yongqi Sun, Neil Saintilan, Senlei Hou, Yanjun Tian, Xieraili Aierken, Yongnian Gao*. HAT-Coast: A multi-source hierarchical adaptive thresholding framework for biennial global mainland shoreline mapping from 2000 to 2020. ISPRS Journal of Photogrammetry and Remote Sensing, 2027, 243P1, 1005-1027.
