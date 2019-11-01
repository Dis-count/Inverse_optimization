
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

    %  4个变量
    X  = new IloNumVar[4];
    for(int i = 0; i < X.length; i++)
      X[i] = cplex.numVar(0.0, 1, IloNumVarType.Int, "X" + i);

    %  初始目标函数
    IloLinearNumExpr obj = cplex.linearNumExpr();
    obj.addTerm(16-8*mu, X[0]);
    obj.addTerm(10-2*mu, X[1]);
    obj.addTerm(0-mu, X[2]);
    obj.addTerm(4-4*mu, X[3]);
    cplex.addMaximize(obj);

    %  约束条件
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
    %  目标函数
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
