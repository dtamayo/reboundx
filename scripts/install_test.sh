#!/bin/bash

# make a new virtual environment (or conda environment) no need to activate
# conda create -n venv pip
# virtualenv venv
# call this with env name as argument, e.g., ./install_test.sh venv

# ***Setup***

# try both virtualenv and conda
source $1/bin/activate
source activate $1

cd ~/
mkdir install_test_temp
cd install_test_temp
git clone https://github.com/hannorein/rebound.git
git clone https://github.com/dtamayo/reboundx.git
cd reboundx
git checkout installtest
cd ~/install_test_temp

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
            f.write(str(e)+'\n')
    except Exception as e:
        f.write('FAILED\n')
        f.write(str(e)+'\n')
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
pip install --no-cache-dir rebound

cd ~/install_test_temp/reboundx
python setup.py clean --all
cd ~/install_test_temp
pip install -i https://testpypi.python.org/pypi --no-cache-dir reboundx

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
pip install --no-cache-dir rebound

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
pip install -i https://testpypi.python.org/pypi --no-cache-dir reboundx

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


