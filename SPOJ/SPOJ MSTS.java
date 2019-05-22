import java.io.OutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.InputMismatchException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Queue;
import java.util.Comparator;
import java.util.LinkedList;
import java.io.InputStream;

/**
 * Built using CHelper plug-in
 * Actual solution is at the top
 *
 * @author Jeel Vaishnav
 */
public class Main {
    public static void main(String[] args) {
        InputStream inputStream = System.in;
        OutputStream outputStream = System.out;
        InputReader in = new InputReader(inputStream);
        PrintWriter out = new PrintWriter(outputStream);
        CountMinimumSpanningTrees solver = new CountMinimumSpanningTrees();
        solver.solve(1, in, out);
        out.close();
    }

    static class CountMinimumSpanningTrees {
        int[] par;
        int[] size;

        int findSet(int i) {
            if (par[i] == i)
                return i;

            par[i] = findSet(par[i]);
            return par[i];
        }

        void union(int i, int j) {
            i = findSet(i);
            j = findSet(j);

            if (i == j)
                return;

            if (size[i] < size[j]) {
                par[i] = j;
                size[j] += size[i];
            } else {
                par[j] = i;
                size[i] += size[j];
            }
        }

        public void solve(int testNumber, InputReader sc, PrintWriter out) {
            int n = sc.nextInt();
            int m = sc.nextInt();

            Edge e[] = new Edge[m];
            for (int i = 0; i < m; ++i)
                e[i] = new Edge(sc.nextInt() - 1, sc.nextInt() - 1, sc.nextInt());

            Arrays.sort(e, new Comparator<Edge>() {

                public int compare(Edge o1, Edge o2) {
                    return o1.weight - o2.weight;
                }
            });

            par = new int[n];
            size = new int[n];
            for (int i = 0; i < n; ++i) {
                par[i] = i;
                size[i] = 1;
            }

            long ans1 = 1L, ans2 = 1L;
            int ptr = 0;
            while (ptr < m) {
                int l = ptr;
                HashMap<Integer, Integer> map = new HashMap<>();
                while (ptr < m && e[ptr].weight == e[l].weight) {
                    if (findSet(e[ptr].u) != findSet(e[ptr].v)) {
                        int u = findSet(e[ptr].u);
                        int v = findSet(e[ptr].v);

                        map.put(u, 1);
                        map.put(v, 1);
                    }

                    ptr++;
                }
                int r = ptr - 1;

                int cur = 0;
                for (int i : map.keySet()) {
                    map.put(i, cur++);
                }

                ArrayList<Integer> adj[] = new ArrayList[cur];
                for (int i = 0; i < cur; ++i)
                    adj[i] = new ArrayList<>();
                for (int i = l; i <= r; ++i) {
                    if (findSet(e[i].u) != findSet(e[i].v)) {
                        int u = map.get(findSet(e[i].u));
                        int v = map.get(findSet(e[i].v));

                        adj[u].add(v);
                        adj[v].add(u);
                    }
                }

                int vis[] = new int[cur];
                int indMap[] = new int[cur];
                int cptr = 0;
                for (int i = 0; i < cur; ++i) {
                    if (vis[i] == 0) {
                        ArrayList<Integer> curIndices = new ArrayList<>();
                        cptr = 0;
                        Queue<Integer> queue = new LinkedList<>();

                        queue.add(i);
                        vis[i] = 1;
                        indMap[i] = cptr++;
                        curIndices.add(i);
                        while (!queue.isEmpty()) {
                            int curInd = queue.poll();

                            for (int j : adj[curInd]) {
                                if (vis[j] == 0) {
                                    queue.add(j);
                                    vis[j] = 1;
                                    indMap[j] = cptr++;
                                    curIndices.add(j);
                                }
                            }
                        }

                        if (cptr == 1)
                            continue;

                        long mat[][] = new long[cptr - 1][cptr - 1];
                        for (int k = 0; k < cptr - 1; ++k) {
                            mat[k][k] = adj[curIndices.get(k)].size();
                            for (int j : adj[curIndices.get(k)]) {
                                if (indMap[j] < cptr - 1)
                                    mat[k][indMap[j]]--;
                            }
                        }

                        ans1 = ans1 * Determinant.findDeterminant(mat, 3L) % 3L;
                        ans2 = ans2 * Determinant.findDeterminant(mat, 10337L) % 10337L;
                    }
                }

                for (int i = l; i <= r; ++i) {
                    union(e[i].u, e[i].v);
                }
            }

            long fans = -1;
            for (int i = 0; i < 3; ++i) {
                if ((10337L * i + ans2) % 3 == ans1)
                    fans = (10337L * i + ans2) % 31011L;
            }

            out.print(fans);
        }

        class Edge {
            int u;
            int v;
            int weight;

            Edge(int a, int b, int c) {
                u = a;
                v = b;
                weight = c;
            }

        }

    }

    static class NumberTheory {
        public static long fast_pow(long a, long b, long mod) {
            long ans = 1L;
            long val = a;

            while (b > 0) {
                if ((b & 1) == 1)
                    ans = ans * val % mod;

                b >>= 1;
                val = val * val % mod;
            }

            return ans;
        }

    }

    static class Determinant {
        public static long findDeterminant(long[][] matrix, long mod) {
            // mat is a square matrix
            // mod is prime

            int n = matrix.length;

            // convert into % mod version and a new mat
            long mat[][] = new long[n][n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (matrix[i][j] < 0) {
                        long req = (-matrix[i][j] - 1L) / mod + 1L; // pigeonhole for ceiling
                        mat[i][j] = (matrix[i][j] + req * mod) % mod;
                    } else
                        mat[i][j] = matrix[i][j];
                    mat[i][j] %= mod;
                }
            }

            long den = 1L;
            long det = 1L;

            for (int col = 0; col < n; ++col) {
                int swapRow = -1;

                for (int row = col; row < n; ++row) {
                    if (mat[row][col] != 0) {
                        swapRow = row;
                        break;
                    }
                }

                if (swapRow == -1)
                    return 0;

                if (swapRow != col) {
                    for (int selectedCol = 0; selectedCol < n; ++selectedCol) {
                        long temp = mat[swapRow][selectedCol];
                        mat[swapRow][selectedCol] = mat[col][selectedCol];
                        mat[col][selectedCol] = temp;
                    }

                    det *= -1L;
                }

                for (int curRow = col + 1; curRow < n; ++curRow) {
                    den = den * mat[col][col] % mod;

                    for (int curCol = n - 1; curCol >= col; --curCol) {
                        mat[curRow][curCol] = (mat[curRow][curCol] * mat[col][col] % mod - mat[col][curCol] * mat[curRow][col] % mod + mod) % mod;
                    }
                }
            }

            if (det < 0)
                det += mod;

            for (int col = 0; col < n; ++col)
                det = det * mat[col][col] % mod;

            det = det * NumberTheory.fast_pow(den, mod - 2, mod) % mod;
            return det;
        }

    }

    static class InputReader {
        private InputStream stream;
        private byte[] buf = new byte[1024];
        private int curChar;
        private int numChars;
        private InputReader.SpaceCharFilter filter;

        public InputReader(InputStream stream) {
            this.stream = stream;
        }

        public int read() {
            if (numChars == -1)
                throw new InputMismatchException();

            if (curChar >= numChars) {
                curChar = 0;
                try {
                    numChars = stream.read(buf);
                } catch (IOException e) {
                    throw new InputMismatchException();
                }

                if (numChars <= 0)
                    return -1;
            }
            return buf[curChar++];
        }

        public int nextInt() {
            int c = read();

            while (isSpaceChar(c))
                c = read();

            int sgn = 1;

            if (c == '-') {
                sgn = -1;
                c = read();
            }

            int res = 0;
            do {
                if (c < '0' || c > '9')
                    throw new InputMismatchException();
                res *= 10;
                res += c - '0';
                c = read();
            }
            while (!isSpaceChar(c));

            return res * sgn;
        }

        public boolean isSpaceChar(int c) {
            if (filter != null)
                return filter.isSpaceChar(c);
            return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == -1;
        }

        public interface SpaceCharFilter {
            public boolean isSpaceChar(int ch);

        }

    }
}

