# Protein3D
An open-source chemistry system, that provides tools to visualize molecular structures and execute pre-built algorithms.


## Install App
- Clone project
    ```
    git clone https://github.com/LuFent/Test-management-system.git
    ```

-  Enable virtualenv (optional)
    - on Linux:
        ```
        virtualenv venv
        source ./venv/bin/activate
        ```
    - on Windows:
        ```
        virtualenv venv
        .\venv\Scripts\activate
        ```
- Install dependencies:
    ```
    pip install -r requirements.txt
    ```

### Local App Setup
```
python main.py
```


### Local Server Start
```
python -m server.main
```
Server will be started at http://localhost:3000/

### Local build

```
python setup.py build
```
Executable file will appear at `Protein3D/build/exe.win-amd64-3.10`:
