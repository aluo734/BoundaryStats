## Notes

### Version 2.1.0 submission

## Test environments
* macOS R 4.3.0
* Rhub
* win-builder

## R CMD check results
* checking for non-standard things in the check directory
    Found the following files/directories:
    ''NULL''
    Found the following files/directories:

* checking for detritus in the temp directory
    'lastMiKTeXException'

Apparent bugs in Rhub.

### Third Submission

> Unknown, possibly misspelled, fields in DESCRIPTION:
   'Author@R'
>
> --> Authors@R
>
> Please fix and resubmit.

Done

 > Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:10.....> or <arXiv:.....>?

Oh, I understand now! Done!

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility

    Maintainer: 'Amy Luo <aluo@vols.utk.edu>'

    New submission

    Possibly misspelled words in DESCRIPTION:\
     Drapeau (16:450) \
     Fortin (16:442) \
     Jacquez (16:463) \
     workflow (16:23)

Workflow is a word missing from the Rhub dictionary, and the others are names.

* checking for non-standard things in the check directory

    Found the following files/directories: 
     ''NULL''
* checking for detritus in the temp directory

   Found the following files/directories:
     'lastMiKTeXException'

These last two notes appear to be bugs in Rhub.

### Second Submission
> Please rather use the Authors@R field and declare Maintainer, Authors and Contributors with their appropriate roles with person() calls.
>
> e.g. something like:
>>  Authors@R: c(person("Alice", "Developer", role = c("aut", "cre","cph"),
                      email = "alice.developer@some.domain.net"),
               person("Bob", "Dev", role = "aut") )

Done

> Please do not start the description with "This package", "Functions to", package name, title or similar.

Done

 > If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
>
>> authors (year) <doi:...>
>>
>> authors (year) <arXiv:...>
>>
>> authors (year, ISBN:...)
>
> or if those are not available: <https:...>
  with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

There are no methods publications attached to this package

> Please add small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.

Done

> \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please replace \dontrun with \donttest.

Done

> Please unwrap the examples if they are executable in < 5 sec, or replace dontrun{} with \donttest{}.

Done

Examples that include boundary_null_distrib() and overlap_null_distrib() take much longer than 5s.

> You write information messages to the console that cannot be easily suppressed.
>
> It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions) -> R/boundary_null_distrib.R; R/overlap_null_distrib.R

Changed cat() to message()

## Test environments
* macOS R 4.3.0
* Rhub
* win-builder

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility

    Maintainer: 'Amy Luo <aluo@vols.utk.edu>'
  
    New submission
  
    Possibly misspelled words in DESCRIPTION:
    workflow (12:23)
  
This isn't a typo , but it doesn't appear to be in the Rhub dictionary.
  
* checking for non-standard things in the check directory
    Found the following files/directories:
    ''NULL''
    Found the following files/directories:

* checking for detritus in the temp directory
    'lastMiKTeXException'

These last two notes appear to be bugs in Rhub.