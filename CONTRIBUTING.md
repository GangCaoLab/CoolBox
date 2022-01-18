## Enviroment

```bash
$ git clone https://github.com/GangCaoLab/CoolBox.git
$ cd CoolBox
$ conda env create --file environment.yml
$ conda activate coolbox
$ python -m pip install . --no-deps -vv
$ python -m pip install ".[doc]"
$ jupyter nbextension enable --py widgetsnbextension
```

## Test 

```bash
$ export PYTHONPATH=`pwd`
$ pytest tests
```

## Rebuild doc

```bash
$ make html
$ make clean
$ make copy
```
