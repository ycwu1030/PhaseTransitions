<!--
 * @Description  : The ISSUES in the Code need to be solved
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-15 19:47:42
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-22 16:44:46
 -->
- [x] #01 Inconsistent between results on Mac and Linux
  - Actually it is a bug in the code, where I didn't keep the variable for `gsl_vector_view`. Under clang, the complier keep it for me, GNU-gcc didn't
- [x] #02 Nan in intermeida nextPoints
  - The nan is due to the minimization problem. I didn't check the convergence status.
- [ ] #03 Need to Accept extra parameters for potential and its derivatives
- [x] #04 GSL insufficient number of points for interpolation error
  - For number of points larger or equal to 3, use the gsl_interp_steffen (GSL >= 2.1)
  - For number of points equal to 2, use the gsl_interp_linear
  - For number of points equal to 1, no interpolation
- [ ] #05 Two phases have just one common point, how to treat such case in removeRedundantPhases
  - For some situation, this case happens when two phases have tiny overlap in temperature due to the poor precision in determining the local minima around that temperature.