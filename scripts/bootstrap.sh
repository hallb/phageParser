#!/usr/bin/env sh
sudo apt install python3 python3-dev python3-venv postgresql postgresql-server-dev-9.5
python3 -m venv venv
# cp scripts/pip.conf venv/pip.conf

echo 'source venv/bin/activate'

