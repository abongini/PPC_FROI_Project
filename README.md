![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
PPC_FROI_Project
</h1>

<br>


## 🔥 Usage

What is it ? 

This project is designed to assist in solving exams questions for the course "Fonctions et Réactions Organiques I," taught by Professor Jieping Zhu. It specifically addresses questions related to both reversible and irreversible carbonyl reactions.

How It Works? 

Upon receiving input in the form of starting materials (i.e., specific molecules and their associated conditions), the application utilizes a pre-defined database of organic reactions to identify the applicable chemical reaction. It then computes and returns the resulting compound, effectively simulating the expected exam response.

Why ? 

 It simplifies identifying reaction pathways and predicting products, thus enhancing study efficiency and accuracy for students preparing their FROI exam.

How ?

You can import and used our package with the following command:
 
```python
from ppc_froi_project import main

# One line to rule them all
main()
```
In addition, this python package returns the name of the reaction if it is present in our database but not coded. Several input cases are possible, such as: starting molecule - condition, first starting molecule - second, first starting molecule - second condition and ending molecule - condition.
This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 

After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n ppc_froi_project python=3.10 
```

```
conda activate ppc_froi_project
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(ppc_froi_project) $ pip install jupyterlab
```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:abongini/PPC_FROI_Project.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(ppc_froi_project) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



