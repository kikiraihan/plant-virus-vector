# Define the path to the virtual environment
venv_dir="/app/venv3114"

cd /app

# Check if the virtual environment exists
if [ -d "$venv_dir" ]; then
  echo "Using existing virtual environment."
else
  echo "Creating a new virtual environment."
  python -m venv $venv_dir
fi

# Activate the virtual environment
source $venv_dir/bin/activate

sudo apt update
sudo apt-get install -y graphviz graphviz-dev

# install requirement pakai poetry
while IFS= read -r line || [ -n "$line" ]; do
    poetry add "$line"
done < requirements_fastapi.txt

# hapus baris pyrdf2vec/__init__.py
# /app/venv3114/lib/python3.10/site-packages/pyrdf2vec/__init__.py
echo "hapus baris nest_asyncio.apply() pyrdf2vec/__init__.py"
sed -i 's/^.*nest_asyncio.apply()/# &/' /app/venv3114/lib/python3.11/site-packages/pyrdf2vec/__init__.py

# Run the application
uvicorn app_fastapi:app --host 0.0.0.0 --port 8009 --reload --log-level debug | tee -a dmy-api.log

# trick: keep container running
while true; do sleep 1; done
