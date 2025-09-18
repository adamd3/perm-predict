FROM node:18-alpine

ARG BACKEND_URL
ENV BACKEND_URL=$BACKEND_URL

WORKDIR /app

COPY package*.json ./

RUN npm install

COPY src ./src
COPY public ./public
COPY next.config.ts .
COPY tsconfig.json .
COPY postcss.config.mjs .
COPY tailwind.config.ts .
COPY components.json .
COPY eslint.config.mjs .

RUN npm run build

EXPOSE 3000

CMD ["npm", "run", "start"]
