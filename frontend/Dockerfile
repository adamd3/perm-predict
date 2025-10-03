FROM node:18-alpine

ARG BACKEND_URL
ENV BACKEND_URL=$BACKEND_URL

WORKDIR /app

COPY package*.json ./

RUN npm install

COPY . .

RUN npm run build

EXPOSE 3000

CMD ["npm", "run", "start"]
