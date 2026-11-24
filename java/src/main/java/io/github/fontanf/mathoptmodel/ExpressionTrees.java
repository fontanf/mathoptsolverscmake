package io.github.fontanf.mathoptmodel;

import java.nio.charset.StandardCharsets;
import java.util.Objects;

/**
 * Expression trees, concatenated, in MathOptModel's layout: pre-order element
 * arrays (root first, each node followed by its children's subtrees) of
 * operator char codes, constant values, variable indices and numbers of
 * children.
 *
 * <p>Operators: {@code '+'} (n-ary sum), {@code '-'} (binary), {@code '*'}
 * (n-ary product), {@code '/'}, {@code 'n'} (negation), {@code 'e'} (exp),
 * {@code 'l'} (log), {@code 'g'} (log10), {@code 'q'} (sqrt), {@code 's'},
 * {@code 'c'}, {@code 't'}, {@code 'h'} (sin, cos, tan, tanh), {@code 'a'}
 * (abs), {@code 'p'} (power), {@code 'm'} / {@code 'M'} (n-ary min / max),
 * {@code 'v'} (variable), {@code 'd'} (defined variable) and {@code 'k'}
 * (constant).
 *
 * <p>The arrays are referenced, not copied.
 */
public final class ExpressionTrees {

    final byte[] operators;
    final int[] numberOfChildren;
    final int[] variables;
    final double[] values;

    /**
     * @param operators operator char codes, one per node
     * @param numberOfChildren number of children of each node
     * @param variables index of each 'v' / 'd' node's (defined) variable, -1 elsewhere; may be null if there are no such nodes
     * @param values value of each 'k' node, 0 elsewhere; may be null (0 everywhere)
     */
    public ExpressionTrees(byte[] operators, int[] numberOfChildren, int[] variables, double[] values) {
        this.operators = Objects.requireNonNull(operators, "operators");
        this.numberOfChildren = Objects.requireNonNull(numberOfChildren, "numberOfChildren");
        this.variables = variables;
        this.values = values;
        Batches.checkLength(numberOfChildren.length, operators.length, "numberOfChildren", "operators");
        if (variables != null)
            Batches.checkLength(variables.length, operators.length, "variables", "operators");
        if (values != null)
            Batches.checkLength(values.length, operators.length, "values", "operators");
    }

    /** Same, with the operators given as a string of ASCII codes, such as {@code "+evpvk"}. */
    public ExpressionTrees(String operators, int[] numberOfChildren, int[] variables, double[] values) {
        this(operators.getBytes(StandardCharsets.US_ASCII), numberOfChildren, variables, values);
    }

    /** The number of nodes. */
    public int numberOfNodes() {
        return operators.length;
    }
}
