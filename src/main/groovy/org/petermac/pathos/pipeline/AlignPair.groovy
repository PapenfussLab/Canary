package org.petermac.pathos.pipeline

import groovy.transform.EqualsAndHashCode

/**
 * Created for PathOS.
 *
 * Description:
 *
 * Hash Class to key the amplicon and both the reads in their entirety
 * The annotation EqualsAndHashCode allows the class created to be used as a HashMap key
 *
 * User: Ken Doig
 * Date: 11-apr2015
 */

@EqualsAndHashCode
class AlignPair
{
    String  ampliconName
    String  read1
    String  read2
}
