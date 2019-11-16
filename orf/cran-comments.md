## Resubmission
* this is the resubmission addressing the comments of Jelena Saf
* thank you for your time to review the package

## Changes in this Resubmission
### Issue 1
* **CRAN**: The LICENSE file is only needed if you have additional restrictions to
the GPL-3 which you have not? In that case omit the file and its
reference in the DESCRIPTION file.
"GPL-3 | file LICENSE" --> "GPL-3"
* **Gabriel**:  I do not have any additional restrictions, so I removed the LICENSE 
file and changed the DESCRIPTION file accordingly. Thereby this issue should be solved.

### Issue 2
* **CRAN**: You write information messages to the console that cannot be easily
suppressed.
Instead of print()/cat() rather use message()/warning()  or
if(verbose)cat(..) if you really have to write text to the console.
(except for print() and summary() functions)
F.i.: margins.R
* **Gabriel**: I checked the usage of print()/cat() in all functions of the package.
The only occurence of print()/cat() outside of print() and summary() functions
is indeed in the file margins.R, namely functions 'margins_output' and 
'margins_output_latex'. However, these two functions are internal only and are
called only from the exported print() and summary() functions. Thus they get
executed only if user calls print() or summary() explicitly and thus fall under
the above mentioned exceptions. Therefore I think this is not an issue.

### Issue 3
* **CRAN**: dontrun{} should be only used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in dontrun{} adds the comment
("# Not run:") as a warning for the user.
Please unwrap the examples if they are executable in < 5 sec, or create
additionally small toy examples to allow automatic testing.
You could also replace dontrun{} with donttest{}, but it would be
preferable to have automatic checks for functions.
* **Gabriel**: Some of the examples take slightly longer than 5 seconds but
are fully executable. Therefore I replaced dontrun{} commands with donttest{}
for these longer examples. Additionally, I removed all the other dontrun{}
commands from all the other examples which are executable under 5 seconds.
These examples demonstrate all the main functions of the package. Thus all 
these examples will get executed in the automatic checks. I think this issue
should be thereby solved.

### Issue 4
* **CRAN**: Following your .Rd documentation, check_discrete_Y() has the same input
and output. This doesn't seem very helpful to me. Please explain.
* **Gabriel**: The function 'check_discrete_Y()' performs checks of the input data from
user. The input vector 'Y' must be discrete and this function throws an error
if it is not or recodes the values of the vector 'Y' if this is neccessary.
F.i. if the input is 'Y=(1,3,5)', this function recodes it to 'Y=(1,2,3)' and 
gives a message about it. This function is internal only and the .Rd documentation 
is internal too. This is analogous to all functions in checks.R. I hope this 
clarifies the issue.

* There were no further issues.
* Additionally to the first submission I have now updated the package version 
in the DESCRIPTION file from '0.1.0' to '0.1.1' and updated the date as well.
* I documented the update in the NEWS.md file.

## Test environments
* local ubuntu 18.04, R 3.6.1
* local macOS 10.15, R 3.6.1
* local windows 10, R 3.6.1
* win-builder (devel and release)

## R CMD check results
* There were 0 ERRORs and 0 WARNINGs and 1 NOTEs. 

* the 1 NOTE is due to the first submission.
* both names were marked as possibly mis-spelled but these are correct.

## Downstream dependencies
* There are currently no downstream dependencies for this package.
