# pastikan saat ini menggunakan versi python 3.8

# Define the path to the virtual environment
venv_dir="venv3816"
version="3.8"

# Check if the virtual environment exists
if [ -d "$venv_dir" ]; then
  echo "Using existing virtual environment."
else
  echo "Creating a new virtual environment."
  python -m venv $venv_dir
fi

# Activate the virtual environment
source $venv_dir/bin/activate

echo "your python version is"
python --version
echo "\n\n"

pip install --upgrade pip
pip install poetry

# Activate again the venv, to make sure poetry install the package in the venv
source ./$venv_dir/bin/activate

# poetry init by pyproject
poetry install

# delete line nest_asyncio.apply() in pyrdf2vec/__init__.py
echo "hapus baris nest_asyncio.apply() di pyrdf2vec/__init__.py"
sed -i '' '/nest_asyncio.apply()/ s/^/#/' $venv_dir/lib/python$version/site-packages/pyrdf2vec/__init__.py