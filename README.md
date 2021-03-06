# 2D Whisker Parameter Analysis

## Introduction
The code is desiged to automatically extract 2D seal whisker parameters and saved on data under *SealWhisker.mat* file, and the program is setup to automatically detect black or green background.  

### Parameters Saved Under SealWhisker.mat
| Parameters Name    | Explanation    |
| :------------- | :------------- |
| Whisker Length       | Length of Seal whisker       |
|row   |   Row number of whisker |
|col   |Col number of Whisker   |
|SealNum   | Harbor number of Whisker   |
|Side   | Left/right side of whisker (1/0) |
|  D_base | Base diameter   |
|D_tip   |Tip diameter   |
|Ratio_R   | D_base/D_tip   |
|whisker_xx   |All x positions alonge the centerline |
|whisker_yy   |All y positions alonge the centerline |
|Valid    |valid or not valid whisker (1/0), it is automatically set to valid value (1), but user can changed value directly from SealWhisker.mat files after running the program   |
|dis_tip_base   |The stright line length from tip to base   |
|Std_upper_concave_xx   | concave x positions on upper side of whisker    |
|Std_upper_concave_yy    | concave y positions on upper side of whisker  |
|Std_upper_convex _xx   | convex x positions on upper side of whisker    |
|Std_upper_convex_yy    | convex y positions on upper side of whisker  |
|Std_upper_concave_xx   | concave x positions on upper side of whisker    |
|Std_upper_concave_yy    | concave y positions on upper side of whisker  |
|Std_upper_convex _xx   | convex x positions on upper side of whisker    |
|Std_upper_convex_yy    | convex y positions on upper side of whisker  |   |   |

## How to use code
### Running new_whisker_parameters.m
- **Step1** :  user needs to edit the whisker infomation
  - Whiker numbers (Row, col, left/right )
  - Harbor number
- **Step2** : crop the ruler image with 11 bars (users need to double clicks the mouse after they cropped the image)
- **Step3** : crop the whisker which one needs to be extracted the parameters (users need to double clicks the mouse after they cropped the image)
- **Step4** : On croppted whisker image, find the whisker's base and tip position by pressing enter key (two times)
  - First time press enter key, the mouth arrow would changed to cross shape
  - The second time press enter key, to locate base/tip position (base must be the first to locate)

### Optional Steps - if you think the whisker is not valid
- **Option 1** : open *SealWhisker.mat* file, and find the *ifvalid* variable, change last row value from 1 to 0;
- **Option 2** : Delete last row for all vairbales under SealWhisker.mat by runing *delet_row.m* file


## Reference image
The follow images are whisker with extracted concave/convex points



![Hierarchy](img/whisker_peaks.png?raw=true)*<center>Figure 2: whisker and centerlines.png</center>*
![Hierarchy](img/rotated_whiskers.png?raw=true)*<center>Figure 2: rotated_whiskers.png</center>*

<!-- ![alt-text-1](img/whisker_peaks.png "title-1") ![alt-text-2](img/whisker_peaks.png "title-2")![alt-text-3](img/whisker_peaks.png "title-3") -->

