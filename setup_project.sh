mkdir perm-predict
cd perm-predict

git init


##------------------------------------------------------------------------------
## Set up frontend
##------------------------------------------------------------------------------

npx create-next-app@latest frontend --typescript --tailwind --eslint




##------------------------------------------------------------------------------
## Set up backend
##------------------------------------------------------------------------------

mkdir backend
cd backend

python3.11 -m venv venv

source venv/bin/activate

pip install fastapi uvicorn python-multipart rdkit scikit-learn
pip freeze > requirements.txt

mkdir app
touch app/main.py



##------------------------------------------------------------------------------
## Start development
##------------------------------------------------------------------------------

cd frontend
npm run dev
