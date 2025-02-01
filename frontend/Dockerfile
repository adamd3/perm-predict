FROM node:18-alpine

WORKDIR /app

COPY package*.json ./

RUN npm install

COPY app/ ./app
COPY public/ ./public
COPY components/ ./components
COPY lib/ ./lib
COPY next.config.js .
COPY tsconfig.json .

RUN npm run build

CMD ["npm", "start"]