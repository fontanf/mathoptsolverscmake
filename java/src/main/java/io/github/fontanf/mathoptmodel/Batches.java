package io.github.fontanf.mathoptmodel;

import java.nio.charset.StandardCharsets;
import java.util.List;

/** Helpers shared by the batch classes. */
final class Batches {

    private Batches() {
    }

    static void checkLength(int length, int expected, String name, String reference) {
        if (length != expected) {
            throw new IllegalArgumentException(name + " has " + length + " entries, but "
                    + reference + " has " + expected + ".");
        }
    }

    static void checkLength(double[] array, int expected, String name, String reference) {
        if (array != null)
            checkLength(array.length, expected, name, reference);
    }

    static void checkLength(int[] array, int expected, String name, String reference) {
        if (array != null)
            checkLength(array.length, expected, name, reference);
    }

    static void checkLength(byte[] array, int expected, String name, String reference) {
        if (array != null)
            checkLength(array.length, expected, name, reference);
    }

    /** Names, as concatenated UTF-8 bytes and offsets (names.size() + 1 entries). */
    record EncodedNames(byte[] data, long[] offsets) {

        static EncodedNames of(List<String> names) {
            if (names == null)
                return new EncodedNames(null, null);
            byte[][] encoded = new byte[names.size()][];
            long[] offsets = new long[names.size() + 1];
            long total = 0;
            for (int i = 0; i < names.size(); ++i) {
                String name = names.get(i);
                if (name == null)
                    throw new IllegalArgumentException("names[" + i + "] is null.");
                encoded[i] = name.getBytes(StandardCharsets.UTF_8);
                offsets[i] = total;
                total += encoded[i].length;
            }
            offsets[names.size()] = total;
            if (total > Integer.MAX_VALUE)
                throw new IllegalArgumentException("names are too long.");
            byte[] data = new byte[(int) total];
            for (int i = 0; i < names.size(); ++i)
                System.arraycopy(encoded[i], 0, data, (int) offsets[i], encoded[i].length);
            return new EncodedNames(data, offsets);
        }
    }
}
