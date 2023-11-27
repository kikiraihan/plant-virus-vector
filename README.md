# plant-virus-vector
this a backend code for [Vektorpedia](https://vektorpedia.ipb.ac.id), an insect vector search engine. This project is connected to [jeki_repo](https://github.com/adibenc/jeki) as the dockerize bundle with frontend.

# instalation
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

## Required folder ignored from git
- certificate_ssl
- venv3114

# about .env
*di server tidak pakai file .env, tapi pakai env yang ada di docker-compose.yml nya.
