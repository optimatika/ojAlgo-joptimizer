/*
 * Copyright 1997-2020 Optimatika
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.ojalgo.joptimizer;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import org.ojalgo.optimisation.Expression;
import org.ojalgo.optimisation.ExpressionsBasedModel;
import org.ojalgo.optimisation.Optimisation;
import org.ojalgo.optimisation.Variable;
import org.ojalgo.structure.Access1D;
import org.ojalgo.structure.Structure1D.IntIndex;
import org.ojalgo.structure.Structure2D.IntRowColumn;

import com.joptimizer.exception.InfeasibleProblemException;
import com.joptimizer.exception.JOptimizerException;
import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.LPOptimizationRequest;
import com.joptimizer.optimizers.LPPrimalDualMethod;
import com.joptimizer.optimizers.NewtonUnconstrained;
import com.joptimizer.optimizers.OptimizationRequest;
import com.joptimizer.optimizers.OptimizationRequestHandler;
import com.joptimizer.optimizers.OptimizationResponse;

import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public final class SolverJOptimizer implements Optimisation.Solver {

    @FunctionalInterface
    public static interface Configurator {

        void configure(final OptimizationRequest request, final OptimizationRequestHandler handler, final Optimisation.Options options);

    }

    static final class Integration extends ExpressionsBasedModel.Integration<SolverJOptimizer> {

        Integration() {
            super();
        }

        public SolverJOptimizer build(final ExpressionsBasedModel model) {

            final boolean max = model.isMaximisation();

            OptimizationRequest request;
            LPOptimizationRequest reqLP;

            if (model.isAnyExpressionQuadratic()) {
                request = new OptimizationRequest();
                reqLP = null;
            } else {
                reqLP = new LPOptimizationRequest();
                request = reqLP;
            }

            final List<Variable> variables = model.getFreeVariables();
            final int numbVars = variables.size();

            final Set<IntIndex> fixed = model.getFixedVariables();

            final Expression objective = model.objective();

            final List<Expression> equalities = model.constraints().filter(expr -> expr.isEqualityConstraint() && !expr.isAnyQuadraticFactorNonZero())
                    .collect(Collectors.toList());
            final int numbEquals = equalities.size();

            if (numbEquals > 0) {

                final DenseDoubleMatrix2D A = new DenseDoubleMatrix2D(numbEquals, numbVars);
                final DenseDoubleMatrix1D b = new DenseDoubleMatrix1D(numbEquals);

                for (int i = 0; i < numbEquals; i++) {
                    final Expression constr = equalities.get(i).compensate(fixed);

                    for (final IntIndex key : constr.getLinearKeySet()) {
                        final int freeIndex = model.indexOfFreeVariable(key);
                        A.set(i, freeIndex, constr.getAdjustedLinearFactor(key));
                    }

                    b.set(i, constr.getAdjustedLowerLimit());
                }

                request.setA(A);
                request.setB(b);
            }

            if (reqLP != null) {
                // LP

                final DenseDoubleMatrix1D c = new DenseDoubleMatrix1D(numbVars);
                final DenseDoubleMatrix1D lb = new DenseDoubleMatrix1D(numbVars);
                final DenseDoubleMatrix1D ub = new DenseDoubleMatrix1D(numbVars);

                for (int i = 0; i < numbVars; i++) {
                    final double cost = objective.getAdjustedLinearFactor(i);
                    c.setQuick(i, max ? -cost : cost);
                    final Variable var = variables.get(i);
                    final double low = var.getUnadjustedLowerLimit();
                    lb.setQuick(i, Double.isFinite(low) ? low : LPPrimalDualMethod.DEFAULT_MIN_LOWER_BOUND);
                    final double upp = var.getUnadjustedUpperLimit();
                    ub.setQuick(i, Double.isFinite(upp) ? upp : LPPrimalDualMethod.DEFAULT_MAX_UPPER_BOUND);
                }

                final List<Expression> lowIneq = model.constraints().filter(expr -> expr.isLowerConstraint() && !expr.isAnyQuadraticFactorNonZero())
                        .collect(Collectors.toList());
                final int numbLowIneq = lowIneq.size();

                final List<Expression> uppIneq = model.constraints().filter(expr -> expr.isUpperConstraint() && !expr.isAnyQuadraticFactorNonZero())
                        .collect(Collectors.toList());
                final int numbUppIneq = uppIneq.size();

                final DenseDoubleMatrix2D G = new DenseDoubleMatrix2D(numbLowIneq + numbUppIneq, numbVars);
                final DenseDoubleMatrix1D h = new DenseDoubleMatrix1D(numbLowIneq + numbUppIneq);

                for (int i = 0; i < numbLowIneq; i++) {
                    final Expression constr = lowIneq.get(i).compensate(fixed);

                    for (final IntIndex key : constr.getLinearKeySet()) {
                        final int freeIndex = model.indexOfFreeVariable(key);
                        G.set(i, freeIndex, -constr.getAdjustedLinearFactor(key));
                    }

                    h.set(i, -constr.getAdjustedLowerLimit());
                }

                for (int i = 0; i < numbUppIneq; i++) {
                    final Expression constr = uppIneq.get(i).compensate(fixed);

                    for (final IntIndex key : constr.getLinearKeySet()) {
                        final int freeIndex = model.indexOfFreeVariable(key);
                        G.set(numbLowIneq + i, freeIndex, constr.getAdjustedLinearFactor(key));
                    }

                    h.set(numbLowIneq + i, constr.getAdjustedUpperLimit());
                }

                reqLP.setC(c);

                reqLP.setLb(lb);
                reqLP.setUb(ub);

                reqLP.setG(G);
                reqLP.setH(h);

            } else {
                // QP

                Expression compObj = objective.compensate(fixed);
                final ConvexMultivariateRealFunction objFunc = SolverJOptimizer.toObjectiveFunction(compObj, numbVars, model);
                request.setF0(objFunc);

                final List<Expression> lowerConstr = model.constraints().filter(expr -> expr.isLowerConstraint()).collect(Collectors.toList());
                final int numbLowerConstr = lowerConstr.size();

                final List<Expression> upperConstr = model.constraints().filter(expr -> expr.isUpperConstraint()).collect(Collectors.toList());
                final int numbUpperConstr = upperConstr.size();

                if ((numbLowerConstr + numbUpperConstr) > 0) {

                    ConvexMultivariateRealFunction[] constrs = new ConvexMultivariateRealFunction[numbLowerConstr + numbUpperConstr];

                    for (int i = 0; i < numbLowerConstr; i++) {
                        Expression compC = lowerConstr.get(0).compensate(fixed);
                        ConvexMultivariateRealFunction c = SolverJOptimizer.toLowerConstraint(compC, numbVars, model);
                        constrs[i] = c;
                    }
                    for (int i = 0; i < numbUpperConstr; i++) {
                        Expression compC = upperConstr.get(0).compensate(fixed);
                        ConvexMultivariateRealFunction c = SolverJOptimizer.toUpperConstraint(compC, numbVars, model);
                        constrs[numbLowerConstr + i] = c;
                    }

                    request.setFi(constrs);
                }

            }

            return new SolverJOptimizer(request, model.options);
        }

        public boolean isCapable(final ExpressionsBasedModel model) {
            return !model.isAnyVariableInteger();
        }

        @Override
        protected boolean isSolutionMapped() {
            return true;
        }

    }

    public static final SolverJOptimizer.Integration INTEGRATION = new Integration();

    static final Configurator DEFAULT = new Configurator() {

        public void configure(final OptimizationRequest request, final OptimizationRequestHandler handler, final Options options) {

            if (request instanceof LPOptimizationRequest) {
                //    ((LPOptimizationRequest) request).setPresolvingDisabled(true);
            }

        }

    };

    private static ConvexMultivariateRealFunction toFunction(final Expression expression, final int dim, boolean negate, final ExpressionsBasedModel model,
            double constant) {

        DenseDoubleMatrix2D p = null;
        DenseDoubleMatrix1D q = null;
        double r = constant;

        if (expression.isAnyQuadraticFactorNonZero()) {
            p = new DenseDoubleMatrix2D(dim, dim);

            for (final IntRowColumn key : expression.getQuadraticKeySet()) {
                final double adjusted = expression.getAdjustedQuadraticFactor(key);
                double d = negate ? -adjusted : adjusted;
                int rowIndex = model.indexOfFreeVariable(key.row);
                int colIndex = model.indexOfFreeVariable(key.column);
                p.set(rowIndex, colIndex, p.get(rowIndex, colIndex) + d);
                p.set(colIndex, rowIndex, p.get(colIndex, rowIndex) + d);
            }

        }

        if (expression.isAnyLinearFactorNonZero()) {
            q = new DenseDoubleMatrix1D(dim);

            for (final IntIndex key : expression.getLinearKeySet()) {
                final double adjusted = expression.getAdjustedLinearFactor(key);
                q.set(model.indexOfFreeVariable(key.index), negate ? -adjusted : adjusted);
            }
        }

        if (p != null) {
            return new PSDQuadraticMultivariateRealFunction(p, q, r);
        } else {
            return new LinearMultivariateRealFunction(q, r);
        }

    }

    static ConvexMultivariateRealFunction toLowerConstraint(final Expression expression, final int dim, final ExpressionsBasedModel model) {
        return SolverJOptimizer.toFunction(expression, dim, false, model, expression.getAdjustedLowerLimit());
    }

    static ConvexMultivariateRealFunction toObjectiveFunction(final Expression expression, final int dim, final ExpressionsBasedModel model) {
        return SolverJOptimizer.toFunction(expression, dim, model.isMaximisation(), model, 0.0);
    }

    static ConvexMultivariateRealFunction toUpperConstraint(final Expression expression, final int dim, final ExpressionsBasedModel model) {
        return SolverJOptimizer.toFunction(expression, dim, false, model, -expression.getAdjustedUpperLimit());
    }

    private final Optimisation.Options myOptions;
    private final OptimizationRequest myRequest;

    SolverJOptimizer(final OptimizationRequest request, final Optimisation.Options options) {
        super();
        myRequest = request;
        myOptions = options;
    }

    public Result solve(final Result kickStarter) {

        final OptimizationRequestHandler handler = myRequest instanceof LPOptimizationRequest ? new LPPrimalDualMethod() : new NewtonUnconstrained(true);

        DEFAULT.configure(myRequest, handler, myOptions);
        final Optional<Configurator> optional = myOptions.getConfigurator(Configurator.class);
        if (optional.isPresent()) {
            optional.get().configure(myRequest, handler, myOptions);
        }

        Optimisation.State retState = Optimisation.State.UNEXPLORED;
        Access1D<Double> retSolution = null;

        if (kickStarter != null) {
            retState = kickStarter.getState();
            myRequest.setInitialPoint(kickStarter.toRawCopy1D());
        }

        handler.setOptimizationRequest(myRequest);

        try {

            handler.optimize();

            final OptimizationResponse response = handler.getOptimizationResponse();

            retState = Optimisation.State.OPTIMAL;
            retSolution = Access1D.wrap(response.getSolution());

        } catch (final InfeasibleProblemException exception) {

            retState = Optimisation.State.INFEASIBLE;

        } catch (final JOptimizerException exception) {

            retState = Optimisation.State.FAILED;

        } finally {

            if ((retSolution == null) && (kickStarter != null)) {
                retSolution = Access1D.wrap(kickStarter.toRawCopy1D());
            }
        }

        return new Optimisation.Result(retState, retSolution);
    }

}
