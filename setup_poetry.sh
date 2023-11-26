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

./poetry_install_requirements.sh

# Run the application
uvicorn app_fastapi:app --host 0.0.0.0 --port 8009 --reload --log-level debug | tee -a dmy-api.log

# trick: keep container running
while true; do sleep 1; done
