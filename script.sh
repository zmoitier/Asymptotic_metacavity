#!/bin/sh

remove_directory() {
    find . -name "$1" -type d \
        -exec echo "removing {}" \; \
        -exec rm -dr {} +
}

remove_file() {
    find . -name "$1" -type f \
        -exec echo "removing {}" \; \
        -exec rm {} +
}

case "$1" in
-c)
    remove_directory ".ipynb_checkpoints"
    remove_directory "__pycache__"
    ;;
-f)
    python3 -m isort .
    python3 -m black .
    ;;
-u)
    python3 -m pip install --user --upgrade pip
    python3 -m pip install --user --upgrade -r requirements.txt
    ;;
*)
    echo "The choice are:"
    echo "  > [-c] for cleaning the temporary python file;"
    echo "  > [-f] for formatting the code;"
    echo "  > [-u] for updating the python package."
    exit 1
    ;;
esac
