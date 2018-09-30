
# Introduction to LazTech Skill-stack

## Preface (Nothing valuable here)

Welcome to Laztech, a dynamically developing astrophysics group led by **Full Professor Alex Lazarian** in University of Wisconsin - Madison.

The **Velocity Gradient Technique (VGT)** is a prospective technique to __map the magnetic field of the universe__ with less cost and higher accuracy.

Don't worry. No one's gonna test you about astrophysics or VGT, nor should you think that having Physics/CS/Statistics background alone can achieve great things on your own. The thing is this: you might need to understand each of them for a little bit, and learn the rest iteratively and cooperatively with our nice yet highly efficient and demanding team members. You will learn a lot -- not just hard skills, but also how to do science like doing a business.

## 1. Basic physics and coding

### Read THIS FIRST:  [Get Started in LazTech Group](https://www.overleaf.com/9279691msypmpcjnpqp#/33470921/).

Read [Get Started in LazTech Group](https://www.overleaf.com/9279691msypmpcjnpqp#/33470921/). The essential basics of physics is stated here. Make sure you read this. Here are several questions:
​      1. What is VGT? Can you explain it in theory and algorithm?
​      2. Can you run the code without error? Why?
​      3. What's the code talks about? What's the result?
​      4. How long did it take? Can it be improved?



## 2. Software Requirements 

These are the softwares that all members will use in their life. Some are covered in the previous section, some are not. 

So, install / try to use these great softwares to improve your coding life:

1. **Julia (0.6.4, and 1.0.0), Python.** Make sure you understand: Basic Julia syntax. Matrix operation. Install/Call Julia module. Call python module. Some libraries in Julia: `HDF5`, `Stats`.Some libraries in Python: `numpy`, `matlibplot`, `astropy`.  (Some questions will be throw later)
2. **Github**. Also learn to use
3. **Anaconda (prefereably 2), Jupyter Notebook / Jupyter Lab**. I know you can just use commandline, but these are great tools that will accelerate your workflow. Besides, @Mike will post some tutorials using the Jupyter Notebook. (If you have trouble install, you can try use these: [Julia Pro](https://juliacomputing.com/products/juliapro.html), Atom (with uber-juno extension)).
4. Good text editors. **[Sublime Text](https://www.sublimetext.com)** and **[Atom](https://atom.io)** will always help you to do some tedious job, especially the *multi-cursor*, *project-wide search & replace* for symbols. even intergrated Julia environment, etc. Learn them in our latter [Question section #2](#Q2 [Simple]: ) to contemplate the beauty of smooth operation. 
5. **Slack**
6. **Teamviewer** and **VNC Viewer** for remote access to our protal.
7. [**Typora**](https://typora.io/) for [markdown](https://markdown-it.github.io/) documentation. You can even use Typora to write some LaTeX, which makes your life much more visualized and easier.
8. **[Overleaf](https://www.overleaf.com/)**. The formal platform to write and generate our LaTeX-based papers. Note that Typora can't do serious LaTeX. Overleaf allows us to share the paper with group members and edit it altogether. Get used to LaTeX and try it.
9. (Win User) [Cygwin](http://www.cygwin.com/) to use a fake-shell on your Windows computer (Mac is great :). Also install an [Ubuntu kernel on Windows](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0) for better experience, performance and thus life.
10. (Mac User) [iTerm2](https://www.iterm2.com/) is a good terminal only if you know how to configurate it. Here's one of the configuration: [iTerm2 + Oh My Zsh + Solarized color scheme + Meslo powerline font](https://gist.github.com/kevin-smets/8568070)

These are a lot. Take a breath and if you have any question, post it on the [issue page](https://github.com/Tairokiya/laztech-intro/issues).



## Questions

or Assignments ...?  You should spend sometime to review these to do things safe, fast and accurately.

### Q1 [Simple]: Write a Julia program that reads all the HDF5 files from a directory, and print out file information. 

Function prototype looks like roughly: `function openHDFFiles (directoryName)`, and you can call the function like this: `openHDFFiles("path/to/some/directory")`. But you can do whatever you want. 

When you done the function writing, download some HDF5 files from google. You can tell us to send you some, but it's good that you explore some too.

#### YouNeedToKnow

1. Julia internal command to change directory `cd()`
2. How to call python library from Julia: `using PyCall; @pyimport some.python.module`
3. Python library function `glob` to grab all things in a directory
4. HDF5 function to show file `h5.open()`

#### Further Question:

1. (**Exception Handling**) How do you prevent the program from interruption when some operation throws exception (Hint: `try-catch ` statement would help you). This is useful when you write a big, overnight program and you don't want the program to stop, but some times things are not going well. (Also, how will `try-catch` influence the performance?)
2. (**Profiling**) How do you print the time of each HDF5 file being read? (Hint: `@time` directory will help you). How do you time the entire for loop? The entire function? What will the output look like? **What tool** can you use to actually see how many times are spend on each line, how many memory is being allocated? **Visualize them!!!**.
3. (**Multiple Return Value**) Make the function return these information: (1) total opened file numbers, (2) every file's size. You should call the function like this: ` filenumbers, filesize = openHDFFiles("path/to/some/directory")`. How should you get filesizes and store them? (Hint: `push!(array, element)`) How should you write the return statement? (Hint: use tuple:`return (filenumbers, filesize_array)`). 
4. (**Coding Style**) There are some general coindg styles ([Google-Java-style](https://google.github.io/styleguide/javaguide.html), [Julia-style guide](https://docs.julialang.org/en/v0.6.2/manual/style-guide/)). Just look through them a little bit. Don't spend too much time on this because the ultimate goal for us is that **your code should be: readable, easily-understood, easy to change, easy to be called by your member**. That's it.



### Q2 [Simple]: Use Sublime Text multi-cursor to transform the file from format A to B

Use Sublime Text to achieve the following things:

#### Format A:

Text Segment 1:

```julia
var1 = 3
somevar2 = 6
a_differnet_var3 = 9
also_a_different_thing = 12
name_it_so_no_comman_subchar_is_shared = 6
```

Text Segment 2:

```
nothing interesting
interested in you
just different thing
so that you have to use
multi-cursor to do it
```

#### Format B:

```julia
var1 = 3 # var1 : nothing interesting
somevar2 = 6 # somevar2 : interested in you
a_differnet_var3 = 9 # a_differnet_var3 : just different thing
also_a_different_thing = 12 # also_a_different_thing : so that you have to use
name_it_so_no_comman_subchar_is_shared = 6 # name_it_so_no_comman_subchar_is_shared : multi-cursor to do it
```

#### YouNeedToKnow

1. Select all the specific character `=` and generate multicursor in your panel (use `Find all`)
2. Move all cursors to the front/end of the line simultaneously
3. Select a word (from that position) of each line, copy and paste them to the end of the line
4. Copy another file of same number of lines and directly concatenate with everyline in this file (simple copy-paste using multi-cursor)

#### Further Questions:

1. (**Indentation**) What should you do if you want to indent every line exactly 4 spaces away from the start of line? (Hint: Select All and then try tab) How about indent 4 spaces *back* to the start of line?

2. ([Regular Expression](https://regex101.com/)101) How to search specific word in the text? How to search for the "start of line"? How to search the "end of line"? 

3. ([Regular Expression](https://regex101.com/)301) How to generate a regex such that it matches a hyperlink? A string wrapped by `""`? A LaTeX formular wrapped by `$$` and might have lines in between? 

4. (**Insert line**) How about insert a blank line after each line? How do you duplicate each line just underneath its original line? 





### Q3[Still-simple] Use Julia to generate some result

This is to get you familar to the modules and algorithms we are using in Julia. If you've heard the 