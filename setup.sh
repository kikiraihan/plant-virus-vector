# source / create new py env if not exist

#!/bin/bash

# Define the path to the virtual environment
# venv_dir="/app/venv_3816"
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

python -m pip install --upgrade pip

# Install required Python packages
# pip install -r requirements.txt
pip install -r r2.txt
sudo apt install graphviz graphviz-dev
pip install pygraphviz

# Run the application
python3 app.py | tee -a dmy-api.log

# trick: keep container running
while true; do sleep 1; done