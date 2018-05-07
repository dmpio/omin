
<!-- ![logo](https://github.com/dmpio/omin/blob/master/images/omin_lil_horns_logo.png) -->


<p align="center">
  <img src="images/omin_lil_horns_logo_reworked_501px_194px.png" >
</p>

<!--
<center>![logo](images/omin_lil_horns_logo_reworked_501px_194px.png)</center>
-->

# omin - Omics Module Integrating and Normalization
### Using python to squeeze the biomarkers out of your big data.
---    

[![submit bug](https://img.shields.io/badge/project%20issues-submit%20bug-red.svg)](https://github.com/dmpio/omin/issues/new?template=issue_template.md&labels=BUG&title=BUG%20:)
[![feature request](https://img.shields.io/badge/project%20issues-submit%20feature%20request-blue.svg)](https://github.com/dmpio/omin/issues/new?template=feature_request.md&labels=FEATURE%20REQUEST&title=FEATURE%20REQUEST%20:)    

[![GitHub issues](https://img.shields.io/github/issues/dmpio/omin.svg)](https://github.com/dmpio/omin/issues)
[![GitHub forks](https://img.shields.io/github/forks/dmpio/omin.svg)](https://github.com/dmpio/omin/network)
[![GitHub stars](https://img.shields.io/github/stars/dmpio/omin.svg)](https://github.com/dmpio/omin/stargazers)

[![PyPI](https://img.shields.io/pypi/v/omin.svg)](https://pypi.org/project/omin/)
![PyPI - Status](https://img.shields.io/pypi/status/omin.svg)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Django.svg)

---
## Requirements

> Anaconda(or Miniconda) -> Python 3.x

[Download Anaconda3](https://docs.anaconda.com/anaconda/install/)

[Download Miniconda3](https://conda.io/miniconda.html)

If you are running windows we suggest that you install Miniconda3 in `%PROGRAMFILES%\Miniconda3` without adding it to your windows path. This way it should not interfere with other installations of Python.

---

## Installation within the Python 3.x environment (Windows 7/8/10)
1. Make sure your computer has Miniconda3/Anaconda3 installed. You can download the free software from the links listed above.
1. Download and unzip the this [zip file](https://github.com/dmpio/omin/archive/master.zip)
1. Open a commandline windows as administrator (press the windows button then type cmd) or teminal (linux users should how to do this, mac users just google it).
1. If Python 3.x is already available on your system path and you know what you are doing skip to step 4.
1.  In the commandline activate Anaconda/Miniconda by typing the following:
    `"%PROGRAMFILES%\Miniconda3\Scripts\activate.bat"`

    or    

    `"%PROGRAMFILES%\Anaconda3\Scripts\activate.bat"`

1. Your commandline now be prepended with the word `(base)`.    
1. `cd` to the unzipped directory from step 1 like so:    

    `cd %userprofile%\Downloads\omin-master\`

1. Then run the command:    

    `pip install .`

1. This command will install the omin package into: `<Your Python distro>\lib\sitepackages`. Now omin can be run from the command-line or imported into jupyter notebook.

---

## Usage:

### Jupyter Notebook/Lab + cookiecutter

1. [`cookiecutter`]() requires `Git`. You can install `Git` if it is not installed already directions can be found here: [installing Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

1. Install [cookiecutter](https://github.com/audreyr/cookiecutter)   
    1. Installing via Anaconda prompt(Windows7/8/10).    
        1. Click the start menu type `anaconda prompt` and press <kbd>Enter</kbd>
        2. Once the prompt has started run the following command:    
            `pip install cookiecutter -U`

1. Using cookiecutter:
    1. Open a command prompt.
    2. Use the [`cd`](#FIXME: Add link) command to change your current directory to where you would like to create your project.
    3. When at your target directory copy and paste the following command then press <kbd>Enter</kbd>:    
        `cookiecutter https://github.com/dmpio/cookiecutter-omin-jupyter-notebook`
    4. `cookiecutter` will prompt you for your name and email then it will create a custom jupyter notebook `.ipynb` file populating it with the information collected from the prompts.
    5. #FIXME: Add a fancy GIF showing behavior.


#### Jupyter Notebook/Lab extensions

1. Install jupyter


---

## Usage: Command line

- Content TK.

---
## Omin overview

<p align="center"><img src="/images/omin_state_diagram.svg" ></p>

---
## Omin process object overview

<p align="center">
  <img src="/images/omin_state_diagram_process_intstance.svg" >
</p>

> State diagrams generated with FreeMind 1.0.1

---
# Contributing

Want to help? Check the [to do list](TO_DO_LIST.md) kill the bugs, make a PR, and bask in the glory<sup>[&#10044;](#asterisk)</sup>.

> <a name="asterisk"> &#10044;</a>: Amount of glory with vary based on bug squashed.

---
<p align="center"><img src="images/duke_octocat_drawing_v1_.300px_292px.png"></p>

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
