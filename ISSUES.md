<!--
 * @Description  : The ISSUES in the Code need to be solved
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-15 19:47:42
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-16 15:30:05
 -->
- [x] #01 Inconsistent between results on Mac and Linux
    - Actually it is a bug in the code, where I didn't keep the variable for `gsl_vector_view`. Under clang, the complier keep it for me, GNU-gcc didn't
- [x] #02 Nan in intermeida nextPoints
    - The nan is due to the minimization problem. I didn't check the convergence status.
- [ ] #03 Need to Accept extra parameters for potential and its derivatives
- [ ] #04 GSL insufficient number of points for interpolation error
- [ ] #05 Two phases have just one common point, how to treat such case in removeRedundantPhases