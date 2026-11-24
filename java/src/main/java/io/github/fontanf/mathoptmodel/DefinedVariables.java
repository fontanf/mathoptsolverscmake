package io.github.fontanf.mathoptmodel;

import java.util.Objects;

/**
 * A batch of defined variables (common subexpressions), for
 * {@link Model#addDefinedVariables(DefinedVariables)}: one expression tree
 * each, referenced from other trees by 'd' nodes. A defined variable may only
 * reference defined variables with a lower index.
 */
public final class DefinedVariables {

    final int[] starts;
    final ExpressionTrees trees;

    /**
     * @param starts offset of each defined variable's tree in {@code trees} (the first one 0)
     * @param trees the trees, concatenated
     */
    public DefinedVariables(int[] starts, ExpressionTrees trees) {
        this.starts = Objects.requireNonNull(starts, "starts");
        this.trees = Objects.requireNonNull(trees, "trees");
    }

    public int number() {
        return starts.length;
    }
}
