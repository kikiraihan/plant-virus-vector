# plant-virus-vector
This is jupyter notebook collections and backend codes for [Vektorpedia](https://vektorpedia.ipb.ac.id), an insect vector search engine based on knowledge graph analysis based.

This repo is connected to dockerize repo [jeki](https://github.com/adibenc/jeki) as backend service container. This repo can be cloned to jeki/src.

# Server instalation
- create venv3114 in this root folder
```
python3 -m venv venv3114
```
- source the venv
```
source ./venv3114/bin/activate
```
- install poetry
```
pip install poetry
```
- run poetry init to create pyproject.toml file
```
poetry init
```
- run setup_fastapi.sh
```
./setup_fastapi.sh
```

# Local instalation
local instalation can use pyproject.toml
- make sure you are using python 3.8.16
- create a file pyproject.toml, you can use pyproject.toml.backup(local3816) as its contents
- run setup_fastapi_local.sh
- run server
```
uvicorn app_fastapi:app --host 0.0.0.0 --port 8009 --reload --log-level debug;
```

## Required folder in server (ignored from git)
- certificate_ssl
- venv3114

# about .env
*di server tidak pakai file .env, tapi pakai env yang ada di docker-compose.yml nya.
