import java.util.*;
import java.text.*;
import java.io.*;
import java.math.*;


public class IVTree {

    public static void main(String[] args) throws Exception {

        double S = 100;
        double K = 105;
        double r = 0.05;
        double q = 0.02;
        double v = 0.25;
        double t = 3.0/12.0;
        String opType = "call";

        int nodes = 4;
        double dt = t / (double) nodes;

        //double put_price = OP(S, K, r, q, v, t, opType, nodes, dt);
        
        // IMPLIED VOLATILITY

        double dv = 0.001;
        double v0 = 0.11, v1 = 1.00, diff, vega;
        double pA = 0, pB = 0, pC = 0;
        double market_price = 3.00;

        while(true){
            pA = OP(S, K, r, q, v0, t, opType, nodes, dt);
            pB = OP(S, K, r, q, v0+dv, t, opType, nodes, dt);
            pC = OP(S, K, r, q, v0-dv, t, opType, nodes, dt);
            vega = (pB - pC)/(2*dv);
            diff = pA - market_price;
            v1 = v0 - diff/vega;
            if(Math.abs(v1 - v0) < dv){
                break;
            }
            v0 = v1;
        }

        double iv_price = OP(S, K, r, q, v1, t, opType, nodes, dt);
        System.out.println(iv_price);
    }

    public static double OP(double S, double K, double r, double q, double v, double t, String optionType, int nodes, double dt){
        
        int row = 4*nodes + 2, col = nodes + 1;
        int center = row / 2 - 1;
        
        double[][] tree = new double[row][col];

        tree[center][0] = S;
        int fp = 2;
        
        double u = U(v, dt);
        double d = D(u);
        double m = M();

        double du = DU(r, q, v, dt);
        double dd = DD(r, q, v, dt);
        double dm = DM(du, dd);
        
        for(int j = 0; j < col; j++){
            for(int i = 1; i < col - j; i++){
                tree[center + i*fp][i + j] = tree[center + (i-1)*fp][j + i - 1]*d;
                tree[center][i + j] = tree[center][j + i-1]*m;
                tree[center - i*fp][i + j] = tree[center - (i-1)*fp][j + i - 1]*u;
            }
        }

        for(int i = 0; i < row; i++){
            if(i % 2 != 0){
                //System.out.println(tree[i][col - 1]);
                if(optionType.equalsIgnoreCase("call")){
                    tree[i][col - 1] = Math.max(tree[i - 1][col - 1] - K, 0);
                } else {
                    tree[i][col - 1] = Math.max(K - tree[i - 1][col - 1], 0);
                }
                    
            }
        }
        

        int cc = 1;
        int jj = 0;
        int iii = 0;

        while(cc < col){
            iii = 3 + jj;
            while(iii < row - jj - 2){
                tree[iii][col - (cc + 1)] = tree[iii][col - cc]*dm + tree[iii - 2][col - cc]*du + tree[iii + 2][col - cc]*dd;
                iii += 2;
            }
            jj += 2;
            cc += 1;
        }

        double option_price = tree[center + 1][0];
        

        return option_price;
    }

    public static void Viz(double[][] tree){
        for(int i = 0; i < tree.length; i++){
            for(int j = 0; j < tree[0].length; j++){
                System.out.print(Math.round(tree[i][j]));
                System.out.print("\t");
            }
            System.out.println();
        }
    }

    public static double DU(double r, double q, double v, double dt){
        double top = Math.exp((r - q)*dt/2.0) - Math.exp(-v*Math.sqrt(dt/2.0));
        double bottom = Math.exp(v*Math.sqrt(dt/2)) - Math.exp(-v*Math.sqrt(dt/2));
        return Math.pow(top/bottom, 2);
    }

    public static double DD(double r, double q, double v, double dt){
        double top = Math.exp(v*Math.sqrt(dt/2.0)) - Math.exp((r - q)*dt/2.0);
        double bottom = Math.exp(v*Math.sqrt(dt/2)) - Math.exp(-v*Math.sqrt(dt/2));
        return Math.pow(top/bottom, 2);
    }

    public static double DM(double du, double dd){
        return 1.0 - (du + dd);
    }
    

    public static double U(double v, double dt){
        return Math.exp(v*Math.sqrt(2.0*dt));
    }

    public static double D(double u){
        return 1.0/u;
    }

    public static double M(){
        return 1.0;
    }

}
