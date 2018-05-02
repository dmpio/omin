
<!-- ![logo](https://github.com/dmpio/omin/blob/master/images/omin_lil_horns_logo.png) -->


<p align="center">
  <img src="images/omin_lil_horns_logo_reworked_501px_194px.png" >
</p>

<!--
<center>![logo](images/omin_lil_horns_logo_reworked_501px_194px.png)</center>
-->

# Omics Modeling Integrating and Normalization(omin)
Using python to squeeze the biomarkers out of your big data.

## Requirements

> Anaconda(or Miniconda) -> Python 3.x

[Download Anaconda3](https://docs.anaconda.com/anaconda/install/)

[Download Miniconda3](https://conda.io/miniconda.html)

If you are running windows we suggest that you install Miniconda3 in `%PROGRAMFILES%\Miniconda3` without adding it to your windows path. This way it should not interfere with other installations of Python.


---

## Installation within the python 3.x environment.
0) Make sure your computer has Miniconda3/Anaconda3 installed. You can download the free software from the links listed above.
1) Download and unzip the this [zip file](https://github.com/dmpio/omin/archive/master.zip)
2) Open a commandline windows as administrator (press the windows button then type cmd) or teminal (linux users should how to do this, mac users just google it).
2.1) If Python 3.x is already available on your system path and you know what you are doing skip to step 4.
3)  In the commandline activate Anaconda/Miniconda by typing the following:
### WINDOWS
```
"%PROGRAMFILES%\Miniconda3\Scripts\activate.bat"
```
or
```
"%PROGRAMFILES%\Anaconda3\Scripts\activate.bat"
```
3.1) Your commandline now be prepended with the word `(base)`.
4) `cd` to the unzipped directory from step 1 like so:
```
cd %userprofile%\Downloads\omin-master\
```
5) Then run the command:

```
pip install .
```
This command will install the omin package into: `<Your Python distro>\lib\sitepackages`

Now omin can be run from the commandline or imported into jupyter notebook.

## Usage: Jupyter Notebook/Lab

- Check out the [omin cookiecutter](https://github.com/dmpio/cookiecutter-omin-jupyter-notebook).

## Usage: Command line


## Omin State Diagram

<p align="center">
  <img src="/images/omin_state_diagram.svg" >
</p>


## Omin Process Object State Diagram

<p align="center">
  <img src="/images/omin_state_diagram_process_intstance.svg" >
</p>

- State diagrams generated with FreeMind 1.0.1

---
<p align="center">
  <img src="images/duke_octocat_drawing_v1_.300px_292px.png">
</p>

## LICENSE:
Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
