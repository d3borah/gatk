package org.broadinstitute.hellbender.metrics;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.Histogram;

import java.util.TreeMap;

/**
 * Custom serializer for htsjdk serializers to work around a kryo bug where classes
 * derived from TreeMap (like htsjdk Histogram) are deserialized and rehydrated as
 * TreeMap instances rather than instances of the derived class.
 */
@SuppressWarnings("rawtypes")
public class HistogramSerializer extends Serializer {

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})
    public void write(Kryo kryo, Output output, Object object) {
        Serializer treeSer = kryo.getDefaultSerializer(TreeMap.class);
        treeSer.write(kryo, output, (TreeMap<?, ?>) object);
        Histogram hist = (Histogram) object;
        output.writeString(hist.getBinLabel());
        output.writeString(hist.getValueLabel());
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})
    public Object read(Kryo kryo, Input input, Class type) {
        Serializer treeSer = kryo.getDefaultSerializer(TreeMap.class);
        TreeMap treeMap = (TreeMap) treeSer.read(kryo, input, TreeMap.class);
        String binLabel;
        String valueLabel;
        binLabel = input.readString();
        valueLabel = input.readString();
        Histogram hist = new Histogram<>(binLabel, valueLabel);
        hist.putAll(treeMap);
        return hist;
    }
}
