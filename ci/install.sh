# Installing miniconda
if [[ "$TRAVIS_OS_NAME" != "windows" ]]; then
    echo "installing miniconda for posix";
    bash $HOME/download/miniconda.sh -b -u -p $MINICONDA_PATH;
elif  [[ "$TRAVIS_OS_NAME" == "windows" ]]; then
    echo "folder $MINICONDA_SUB_PATH does not exist"
    echo "installing miniconda for windows";
    choco install miniconda3 --params="'/JustMe /AddToPath:1 /D:$MINICONDA_PATH_WIN'";
fi;

export PATH="$MINICONDA_PATH:$MINICONDA_SUB_PATH:$MINICONDA_LIB_BIN_PATH:$PATH";

# Checking miniconda existance
echo "checking if folder $MINICONDA_SUB_PATH exists"
if [[ -d $MINICONDA_SUB_PATH ]]; then
    echo "folder $MINICONDA_SUB_PATH exists"
else
    echo "folder $MINICONDA_SUB_PATH does not exist"
fi;

source $MINICONDA_PATH/etc/profile.d/conda.sh;
hash -r;
echo $TRAVIS_OS_NAME
echo $CONDA_PYTHON
python --version
conda config --set always_yes yes --set changeps1 no;
conda update -q conda;

# Free channel has Python < 3.7
conda config --set restore_free_channel true;

# Useful for debugging any issues with conda
conda info -a

# See if test-environment already available
echo "create test-environment";
conda create -n test-environment;

conda activate test-environment
conda list
