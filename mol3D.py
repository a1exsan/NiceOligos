from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from nicegui import ui
from rdkit.Chem import rdDistGeom

# Генерация молекулы
smiles = """
Cc1cn([C@H]2C[C@H](OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(N)nc4=O)C[C@@H]3O)[C@@H](COP(=O)(O)O[C@H]3C[C@H](n4ccc(N)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4ccc(N)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cnc5c(N)ncnc54)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4ccc(N)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4ccc(N)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cc(C)c(O)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cc(C)c(O)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cnc5c(N)ncnc54)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cc(C)c(O)nc4=O)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cnc5c(N)ncnc54)O[C@@H]3COP(=O)(O)O[C@H]3C[C@H](n4cnc5c(N)ncnc54)O[C@@H]3COP(=O)(O)OC3CCC(NC(=O)CCCc4cn(CCC/N=C(/O)c5cc(Cl)c6c(c5Cl)C5(OC6=O)c6cc(-c7ccccc7)c(O)cc6Oc6cc(O)c(-c7ccccc7)cc65)nn4)CC3)O2)c(=O)nc1O
"""
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

params = rdDistGeom.EmbedParameters()
params.clearConfs = True         # пример параметров
params.randomSeed = 42
params.maxIterations = 1000
#params.pruneRmsThresh = 0.1
conf_id = AllChem.EmbedMolecule(mol, params)

#res = AllChem.EmbedMolecule(mol, maxAttempts=1000, randomSeed=42)
print(conf_id)
if conf_id == -1:
    print("Генерация конформеров не удалась")
else:
    AllChem.UFFOptimizeMolecule(mol)
mol_block = Chem.MolToMolBlock(mol)

# Создаем py3Dmol view
view = py3Dmol.view(width=1500, height=1000)
view.addModel(mol_block, 'mol')
view.setStyle({'stick': {}})
view.zoomTo()
html = view._make_html()

# Разделяем скрипты и контейнер
# В py3Dmol html полный html-документ, значит, выделим body содержимое и добавим через add_body_html

# Здесь для примера добавим весь созданный html как is_body_html
ui.add_body_html(html)

# Добавляем только div-контейнер для визуализации (ID должен совпадать с id py3Dmol)
ui.html('<div id="container" style="width: 1500px; height: 1000px;"></div>').style('display: block;')
ui.run(port=8085)