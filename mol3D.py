from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from nicegui import ui

# Генерация молекулы
smiles = """
COc1cc(/N=N/c2ccc(C)cc2[N+](=O)[O-])c(C)cc1/N=N/c1ccc(N(CCO)CCOP(=O)(O)O[C@H]2C[C@H](n3ccc(N)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(O)nc(N)nc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(O)nc(N)nc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3ccc(N)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(O)nc(N)nc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3ccc(N)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cc(C)c(O)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(O)nc(N)nc43)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3ccc(N)nc3=O)O[C@@H]2COP(=O)(O)O[C@H]2C[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)OCCCCCC/N=C(\O)c2ccc3c(c2)C2(OC3=O)c3ccc(O)cc3Oc3cc(O)ccc32)cc1
"""
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
#AllChem.UFFOptimizeMolecule(mol)
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