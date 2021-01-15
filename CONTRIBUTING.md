## Enviroment

``bash
git clone https://github.com/GangCaoLab/CoolBox.git
cd CoolBox
conda env create --file environment.yml
conda activate coolbox
python -m pip install . --no-deps -vv
jupyter nbextension enable --py widgetsnbextension
```

## Test and rebuild doc

```bash
pytest tests
cd docs && make html; cp -r build/html/* .
```
