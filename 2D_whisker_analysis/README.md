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
- **Step1** :  user need to edit the whisker infomation on the next section
  - Whiker numbers (Row, col, left/right )
  - harbor number
- **Step2** : crop the ruler image with 11 bars (use needs to double clicks after they cropped Image)
- **Step3** : crop the whisker which one needs to be extracted the parameters (user need to double clicks after they cropped Image)
- **Step4** : On croppted whisker image, find the whisker base and tip position by pressing enter key (two times)
  - First time press enter key, the mouth arrow would changed to crop shape
  - The second time press enter key, to locate base/tip position (base must be the first to locate)

### Optional Steps - if you think the whisker is not valid
- **Option 1** : open SealWhisker.mat file, and find ifvalid variable, change last row from 1 to 0;
- **Option 2** : Delete last row for all vairbales under SealWhisker.mat by *runing delet_row.m* file


## Reference image
The follow image is a whisker with extracted concave/convex points

![whisker_peaks.png](https://github.com/SeNSE-lab/robots_sealwhiskers/blob/master/2D_whisker_analysis/img/whisker_peaks.png)![rotated_whiskers.png](https://github.com/SeNSE-lab/robots_sealwhiskers/blob/master/2D_whisker_analysis/img/rotated_whiskers.png)
