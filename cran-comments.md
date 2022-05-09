# Resubmission
This is a resubmission. In this version I have:

* Replaced all "print()" function calls with "message()"

* Reset user's "par()" after changing them in example code.
  
## R CMD check results
### Platform:	Ubuntu Linux 20.04.1 LTS, R-release, GCC
There were no ERRORs or WARNINGs. There was one note.
> Possibly misspelled words in DESCRIPTION:
  > Baldanzini (10:157, 10:523)  
  
  last name of an author.
  
  > Pierini (10:173, 10:539)
  
  last name of an author.
  
  > SQG (10:101, 10:679)
  
  abbreviation for an algorithm.
    
  > Sala (10:151, 10:517) 
  
  last name of an author.
    
  > optim (10:17)
  
  name of a function in base R.

### Platform:	Fedora Linux, R-devel, clang, gfortran
There were no ERRORs or WARNINGs. There was one note.
Same as above

### Platform:	Windows Server 2022, R-devel, 64 bit
There were no ERRORs or WARNINGs. There was two note.
Same as above + 
> * checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

The second note might be due to a bug/crash in MiKTeX, as is noted in R-hub issue #503, which can be found at https://github.com/r-hub/rhub/issues/503.
