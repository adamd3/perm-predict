FROM node:18-alpine

ARG BACKEND_URL
ENV BACKEND_URL=$BACKEND_URL

WORKDIR /app

COPY package*.json ./

RUN npm install

# Create the public directory
RUN mkdir -p ./public

# Move RDKit WASM and JS files from node_modules to the public directory
RUN cp node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm ./public/RDKit_minimal.wasm
RUN cp node_modules/@rdkit/rdkit/dist/RDKit_minimal.js ./public/RDKit_minimal.js

COPY . .

RUN npm run build

EXPOSE 3000

CMD ["npm", "run", "start"]
