# Development Environment Setup Scripts

1. ```scripts/bootstrap.sh``` uses apt to get required software including python, postgres, etc...
2. ```scripts/server.sh``` starts the postgres service
3. ```source venv/bin/activate``` to start virtualenv
4. ```scripts/setup.sh``` uses pip to get python dependencies, creates and migrates database structure

