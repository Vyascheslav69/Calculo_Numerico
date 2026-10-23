import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

public class RegressaoQuadratica {

    private static final int MINIMO_DE_PONTOS = 4;

    public static void main(String[] args) {
        String arquivo = (args.length > 0) ? args[0] : "dados.csv";

        try {
            // ler os pontos do arquivo CSV.
            double[][] pontos = lerCsv(arquivo);
            double[] x = pontos[0];
            double[] y = pontos[1];

            // montar e resolver as equações normais.
            double[] alfas = calcularRegressaoQuadratica(x, y);

            // imprimir α2, α1 e α0.
            imprimirSolucao(alfas);

        } catch (IOException e) {
            System.err.println("Erro ao ler o arquivo " + arquivo + ": " + e.getMessage());
            System.exit(1);
        } catch (IllegalArgumentException e) {
            System.err.println("Erro: " + e.getMessage());
            System.exit(1);
        }
    }

    // LEITURA DO CSV
    static double[][] lerCsv(String caminho) throws IOException {
        List<double[]> pontos = new ArrayList<>();

        try (BufferedReader leitor = new BufferedReader(new FileReader(caminho))) {
            String linha;
            int numeroDaLinha = 0;

            while ((linha = leitor.readLine()) != null) {
                numeroDaLinha++;
                linha = linha.trim();

                if (linha.isEmpty()) {
                    continue;
                }

                // Separa o x do f(x) usando a vírgula.
                String[] partes = linha.split(",");
                if (partes.length != 2) {
                    throw new IllegalArgumentException("linha " + numeroDaLinha
                            + " deveria ter o formato x,f(x): \"" + linha + "\"");
                }

                try {
                    double xi = Double.parseDouble(partes[0].trim());
                    double yi = Double.parseDouble(partes[1].trim());
                    pontos.add(new double[] {xi, yi});
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("linha " + numeroDaLinha
                            + " contém um valor que não é número: \"" + linha + "\"");
                }
            }
        }

        // Copia a lista para dois vetores simples (x e y)
        int n = pontos.size();
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = pontos.get(i)[0];
            y[i] = pontos.get(i)[1];
        }
        return new double[][] {x, y};
    }

    // CÁLCULO DOS COEFICIENTES α2, α1, α0
    /* Calcula α2, α1 e α0 da parábola φ(x) = α2·x² + α1·x + α0 que melhor se ajusta aos pontos
     * Retorna o vetor {α2, α1, α0}. */

    static double[] calcularRegressaoQuadratica(double[] x, double[] y) {
        int n = x.length;

        // O trabalho exige no mínimo 4 pontos.
        if (n < MINIMO_DE_PONTOS) {
            throw new IllegalArgumentException("são necessários no mínimo "
                    + MINIMO_DE_PONTOS + " pontos, mas o arquivo tem " + n + ".");
        }

        // Uma parábola só fica bem definida se existirem pelo menos 3 valores de x DIFERENTES.
        Set<Double> xDistintos = new HashSet<>();
        for (double xi : x) {
            xDistintos.add(xi);
        }
        if (xDistintos.size() < 3) {
            throw new IllegalArgumentException("são necessários pelo menos 3 valores"
                    + " de x diferentes para definir uma parábola única.");
        }

        // Calcular os somatórios
        // Percorre os pontos uma vez e acumula todas as somas que aparecem nas equações normais.
        double somaX = 0, somaX2 = 0, somaX3 = 0, somaX4 = 0;
        double somaY = 0, somaXY = 0, somaX2Y = 0;

        for (int i = 0; i < n; i++) {
            double xi = x[i];
            double yi = y[i];

            // Reaproveitar as potências já calculadas
            double xi2 = xi * xi;     // x²
            double xi3 = xi2 * xi;    // x³ = x² · x
            double xi4 = xi2 * xi2;   // x⁴ = x² · x²

            somaX   += xi;
            somaX2  += xi2;
            somaX3  += xi3;
            somaX4  += xi4;

            somaY   += yi;
            somaXY  += xi * yi;
            somaX2Y += xi2 * yi;
        }

        // montar a matriz M e o vetor B (M · α = B)
        double[][] M = {
            {somaX4, somaX3, somaX2},
            {somaX3, somaX2, somaX },
            {somaX2, somaX,  n     }
        };
        double[] B = {somaX2Y, somaXY, somaY};

        // PASSO 4: resolver o sistema pela Regra de Cramer
        return resolverPorCramer(M, B);
    }

    // Resolve um sistema linear 3x3  M · α = B  pela REGRA DE CRAMER.
    static double[] resolverPorCramer(double[][] M, double[] B) {
        double detM = determinante3x3(M);

        // Se det(M) = 0 não existe solução única
        if (detM == 0.0) {
            throw new IllegalArgumentException("o sistema não possui solução única"
                    + " (determinante igual a zero).");
        }

        double[] solucao = new double[3];
        for (int coluna = 0; coluna < 3; coluna++) {
            double[][] Mj = trocarColuna(M, coluna, B);
            solucao[coluna] = determinante3x3(Mj) / detM;
        }
        return solucao; // {α2, α1, α0}
    }

    // Devolve cópia da matriz M com a coluna indicada substituída pelo vetor B.
    private static double[][] trocarColuna(double[][] M, int coluna, double[] B) {
        double[][] copia = new double[3][3];
        for (int linha = 0; linha < 3; linha++) {
            copia[linha] = M[linha].clone();
            copia[linha][coluna] = B[linha];
        }
        return copia;
    }

    // Determinante de uma matriz 3x3 pela expansão de Laplace
    private static double determinante3x3(double[][] m) {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // print solução
    static void imprimirSolucao(double[] alfas) {
        System.out.printf(Locale.US, "%.6f%n", alfas[0]); // α2
        System.out.printf(Locale.US, "%.6f%n", alfas[1]); // α1
        System.out.printf(Locale.US, "%.6f%n", alfas[2]); // α0
    }
}
