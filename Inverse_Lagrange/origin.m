
package lagranger;

import java.io.IOException;
import ilog.concert.IloException;

public class MainFrame {
  double best_ub;
  double best_lb;
  double best_mu;
  double[] best_sl;
  Subproblem sp;
  public MainFrame()
  {
    best_lb = 0;
    best_ub = 1e10;
    sp = new Subproblem();
    best_sl = new double[4];
  }

  // 次梯度方法求解拉格朗日对偶
  public void solve(double min_step_size, int max_iter) throws IOException, IloException
  {
    int iter = 0;
    int non_improve = 0;
    int max_non_improve = 3;
    double lambda = 2;
    double step_size = 1;
    double mu = 0;         // 初始化拉格朗日乘子
    sp.construct(mu);      // 松弛第一个约束条件的拉格朗日松弛
    while(iter++ < max_iter)
    {
      sp.changeObj(mu);
      if (sp.solve() == false)
      {
        System.out.println("The Lagrangian problem solve wrong!");
        System.exit(0);
      }

      // 更新上界
      if(sp.opt_cost < best_ub)
      {
        best_ub = sp.opt_cost;
        best_mu = mu;
        for(int i = 0; i < best_sl.length; i++)
          best_sl[i] = sp.opt_x[i];
        non_improve = 0;
      }
      else
        non_improve++;
      System.out.println("iter " + iter + "******************************");
      System.out.println("best lb " + best_lb);
      System.out.println("best ub " + best_ub);
      System.out.println("current ub " + sp.opt_cost);
      System.out.println("mu " + mu);
      double subgradient = 8*sp.opt_x[0] + 2*sp.opt_x[1] + sp.opt_x[2] + 4*sp.opt_x[3] - 10;

      mu = Math.max(0, mu + step_size * subgradient);

      // 满足原问题约束的可行解可以作为原问题的下界
      if (subgradient <= 0)
      {
        double current_lb = 16*sp.opt_x[0] + 10*sp.opt_x[1] + 4*sp.opt_x[3];
        if (current_lb > best_lb)
          best_lb = current_lb;
      }

      // 上界未更新达到一定次数
      if(non_improve >= max_non_improve)
      {
        lambda /= 2;
        non_improve = 0;
      }

      double dist = Math.pow(subgradient, 2);

      // 迭代停止条件2和3
      if(dist <= 0.0 || best_lb >= best_ub - 0.0000001)
        break;

      step_size = lambda * (sp.opt_cost - best_lb) / dist;

      // 迭代停止条件4
      if(step_size < min_step_size)
        break;
    }
  }
  public static void main(String[] args) throws IOException, IloException
  {
    MainFrame mf = new MainFrame();
    mf.solve(0.01, 10);
    System.out.println("result: ");
    System.out.println("best_lb: " + mf.best_lb);
    System.out.println("best_ub: " + mf.best_ub);
    double gap = Math.round((mf.best_ub - mf.best_lb) *  10000 / mf.best_ub) / 100;
    System.out.println("gap: " + gap + "%");
  }
}






package lagranger;

import ilog.concert.*;
import ilog.cplex.IloCplex;

public class Subproblem {
  IloCplex cplex;
  double opt_cost;
  double mu;
  double[] opt_x;
  IloNumVar[] X;
  public void construct(double cmu) throws IloException
{
    cplex = new IloCplex();
    cplex.setOut(null);
    mu = cmu;

    // 4个变量
    X  = new IloNumVar[4];
    for(int i = 0; i < X.length; i++)
      X[i] = cplex.numVar(0.0, 1, IloNumVarType.Int, "X" + i);

    // 初始目标函数
    IloLinearNumExpr obj = cplex.linearNumExpr();
    obj.addTerm(16-8*mu, X[0]);
    obj.addTerm(10-2*mu, X[1]);
    obj.addTerm(0-mu, X[2]);
    obj.addTerm(4-4*mu, X[3]);
    cplex.addMaximize(obj);

    // 约束条件
    IloLinearNumExpr expr1 = cplex.linearNumExpr();
    expr1.addTerm(1, X[0]);
    expr1.addTerm(1, X[1]);
    cplex.addLe(expr1, 1);
    IloLinearNumExpr expr2 = cplex.linearNumExpr();
    expr1.addTerm(1, X[2]);
    expr1.addTerm(1, X[3]);
    cplex.addLe(expr2, 1);
  }

  public void changeObj(double cmu) throws IloException
{
    // 目标函数
    mu = cmu;
    IloLinearNumExpr new_obj = cplex.linearNumExpr();
    new_obj.addTerm(16-8*mu, X[0]);
    new_obj.addTerm(10-2*mu, X[1]);
    new_obj.addTerm(0-mu, X[2]);
    new_obj.addTerm(4-4*mu, X[3]);
    cplex.getObjective().clearExpr();
    cplex.getObjective().setExpr(new_obj);
  }

  public boolean solve() throws IloException
{
    if(this.cplex.solve())
    {
      opt_cost = cplex.getObjValue() + 10*mu;
      opt_x = new double[X.length];
      for (int i = 0; i < X.length; i++)
        opt_x[i] = cplex.getValue(X[i]);
      return true;
    }
    cplex.exportModel("model.lp");
    return false;
  }
}
