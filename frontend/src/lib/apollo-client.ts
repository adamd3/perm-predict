import { ApolloClient, InMemoryCache, createHttpLink } from '@apollo/client';

const httpLink = createHttpLink({
  uri: process.env.NEXT_PUBLIC_GRAPHQL_ENDPOINT || '/api/graphql',
});

console.log("NEXT_PUBLIC_GRAPHQL_ENDPOINT:", process.env.NEXT_PUBLIC_GRAPHQL_ENDPOINT);

export const apolloClient = new ApolloClient({
  link: httpLink,
  cache: new InMemoryCache(),
  defaultOptions: {
    watchQuery: {
      pollInterval: 2000, // Poll every 2 seconds for job results
    },
  },
});