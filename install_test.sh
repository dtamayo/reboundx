#!/bin/bash

# make and activate a new virtual environment (or conda environment) with pip installed
# conda create -n venv pip
# virtualenv venv
# set venv to the name of the environment

venv=t4

# ***Setup***

# try both virtualenv and conda
source $venv/bin/activate
source activate $venv

cd ~/
mkdir install_test_temp
cd install_test_temp
git clone https://github.com/dtamayo/rebound.git
git clone https://github.com/dtamayo/reboundx.git
cd rebound
git checkout pypi
cd ../reboundx
git checkout pypi

pip uninstall -y reboundx
pip uninstall -y rebound

# ***Script to run each time***
pyscript="
with open('results.txt', 'a') as f:
	try:
		import reboundx
		exception = reboundx.install_test()
		if exception is None:
			f.write('OK\n')
		else:
			f.write('FAILED\n')
			f.write(e)
	except Exception as e:
		f.write('FAILED\n')
		f.write(str(e))
"

# ***Run all combinations***

setup="REBOUND: Local REBOUNDx: Local: "

cd ~/install_test_temp/rebound
python setup.py clean --all
pip install -e .

cd ~/install_test_temp/reboundx
python setup.py clean --all
pip install -e .

cd ~/install_test_temp
python -c "
with open ('results.txt', 'a') as f:
	f.write('$setup')
"
python -c "$pyscript"

pip uninstall -y reboundx
pip uninstall -y rebound

setup="REBOUND: PyPI REBOUNDx: PyPI: "

cd ~/install_test_temp/rebound
python setup.py clean --all
cd ~/install_test_temp
pip install -i https://testpypi.python.org/pypi rebound

cd ~/install_test_temp/reboundx
python setup.py clean --all
cd ~/install_test_temp
pip install -i https://testpypi.python.org/pypi reboundx

cd ~/install_test_temp
python -c "
with open ('results.txt', 'a') as f:
	f.write('$setup')
"
python -c "$pyscript"

pip uninstall -y reboundx
pip uninstall -y rebound

setup="REBOUND: PyPI REBOUNDx: Local: "

cd ~/install_test_temp/rebound
python setup.py clean --all
cd ~/install_test_temp
pip install -i https://testpypi.python.org/pypi rebound

cd ~/install_test_temp/reboundx
python setup.py clean --all
pip install -e .

cd ~/install_test_temp
python -c "
with open ('results.txt', 'a') as f:
	f.write('$setup')
"
python -c "$pyscript"

pip uninstall -y reboundx
pip uninstall -y rebound

setup="REBOUND: Local REBOUNDx: PyPI: "

cd ~/install_test_temp/rebound
python setup.py clean --all
pip install -e .

cd ~/install_test_temp/reboundx
python setup.py clean --all
cd ~/install_test_temp
pip install -i https://testpypi.python.org/pypi reboundx

cd ~/install_test_temp
python -c "
with open ('results.txt', 'a') as f:
	f.write('$setup')
"
python -c "$pyscript"

pip uninstall -y reboundx
pip uninstall -y rebound

cat ~/install_test_temp/results.txt

cd ~/
rm -rf install_test_temp


