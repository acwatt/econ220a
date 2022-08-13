Using the following PDF, I needed to make several syntax changes to the old Matlab 5.1 code for it to run on Matlab R2022a.
"Optimization Toolbox -- User�s Guide, Version 2"
https://instruct.uwo.ca/engin-sc/391b/downloads/optim_tb.pdf

Pg 96, section "Converting Your Code to Version 2 Syntax" has a syntax convertion table that I used to make the following changes to the code:

# fminu -> fminunc
rc_dc.m, line 78 needs to be updated.
```
[theta2, options] = fminu('gmmobj',theta2,options,'gradobj')
```

Pg 101 says:
> In Version 1.5, you used this call to fminu.
```
[X,OPTIONS] = fminu('FUN',x0,OPTIONS,'GRADFUN',P1,P2,...);
```
> with `F = FUN(X,P1, ...)` and `G = GRADFUN(X,P1, ...)`.

> ...

> If you have an existing FUN and GRADFUN that you do not want to rewrite, you can pass them both to fminunc by placing them in a cell array.
```
OPTIONS = optimset('GradObj','on'); % Gradient is supplied
[X,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc({@FUN,@GRADFUN},x0,OPTIONS,P1,P2,...);
```

Comparing to our line 78, we need to change `fminu -> fminunc` from our code. Also note that the example variable/function names map to our names as:
- `X -> theta2`
- `OPTIONS -> options`
- `FUN -> gmmobj`
- `x0 -> theta2`
- `GRADFUN -> gradobj`


So I changed
```
[theta2, options] = fmin('gmmobj',theta2,options,'gradobj')
```
to
```
options = optimset('GradObj','on'); % Gradient is supplied
[theta2,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc({@gmmobj,@gradobj},theta2,options);
```
Note that 'GradObj' is just a string used by `optimset` to create an options object. This does not need to be changed. Also note that `theta2` is both an input and output -- it is being updated from the initial values to more optimal values. The above two lines replaced line 78.





# option(8) -> FVAL
Because of the change above to `fminu`, the variable `options` no longer contains the minimized GMM function value. It is now stored in `FVAL`. So I updated the display lines at the end of rc_dc.m (lines 118-122). Updated `num2str(option(8))` to `num2str(FVAL)`





# option(10) -> OUTPUT.iterations
Same as above, updated `num2str(option(10))` to `num2str(OUTPUT.iterations)`





# semcoef sparse matric -> full matrix
For displaying at the end only, I had to use `full(semcoef(i))` instead of just `semcoef(i)` because `semcoef` is a sparce matrix and was displaying not in the intended format. `full()` forces it to display like a usual matrix.




Any modified files are subscripted with "_acw". All other files remain in their original state.