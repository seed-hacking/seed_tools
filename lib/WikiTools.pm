#!/usr/bin/perl -w

#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

package WikiTools;

    use strict;
    use Tracer;
    use LWP::UserAgent;
    use FIG_Config;

=head1 Wiki-Handling Methods

=head2 Introduction

This package handles functions related to the NMPDR Wiki. It is also
used by [[Erdbpm]] to generate Wiki output. It can be replaced by a different
package if documentation output of a different type is desired.

The fields in this object are as follows.

=over 4

=item userAgent

An LWP object that can be used to send requests to the Wiki server.

=item username

User name for logging on to the wiki server.

=item password

Password for logging on to the wiki server.

=item url

URL of the wiki server access script.

=back

=cut

=head2 Public Methods

=head3 new

   my $wiki = WikiTools->new();

Construct a new Wiki object.

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Create the Wiki-Handling Methods object.
    my $retVal = {
                    username => $FIG_Config::nmpdr_wiki->{username},
                    password => $FIG_Config::nmpdr_wiki->{password},
                    url => $FIG_Config::nmpdr_wiki->{url}
                };
    # Add the user agent.
    $retVal->{userAgent} = LWP::UserAgent->new();
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 WebApplicationObject

    my $appObject = WikiTools::WebApplicationObject($session);

This method returns a [[WebApplication]] application object for the Wiki.
It is used by the authentication and access control methods to access the
[[WebAppBackend]] database.

=cut

sub WebApplicationObject {
    # Get the parameters.
    my ($session) = @_;
    # Declare the return variable.
    my $retVal = DBMaster->new(-database => $FIG_Config::webapplication_db,
                               -host => $FIG_Config::webapplication_host,
                               -user => $FIG_Config::webapplication_user);
    # Return the result.
    return $retVal;
}

=head3 ComputeTitle

    my ($title, $project) = WikiTools::ComputeTitle($fileName, $path);

This method computes the project name and Wiki title for a source file.
The specified path must be relative to the distribution directory. So,
for example, the path for the Ajax component would be
C<WebApplication/WebComponent>. For almost every source file, the path
will be a single word (e. g. C<FigKernelPackages>). The project name is
the first directory in the incoming path, and the Wiki title is computed
by converting the file name to capital case and squeezing out the
underscores.

=over 4

=item fileName

Unqualified name of the specified source file (e.g. C<FIG.pm>).

=item path

Path in the source file tree to the specified file.

=item RETURN

Returns a two-element list. The first element is a WikiWord for the file
title and the second is a project name.

=back

=cut

sub ComputeTitle {
    # Get the parameters.
    my ($fileName, $path) = @_;
    # Compute the project name from the project path.
    my $project = ($path =~ /^(\w+)/ ? $1 : $path);
    # Turn the project name into a wiki word.
    $project = Wikify($project, "Project");
    # Compute the file title from the file name.
    my $title = Wikify($fileName);
    # Return the results.
    return ($title, $project);
}

=head3 Wikify

    my $wikiword = WikiTools::Wikify($string, $suffix);

Convert a string into a [[WikiWord]]. WikiWords are in [[capital case]],
cannot contain underscores, and must have at least two humps (that is,
two transitions from a capital letter to a digit or small letter). To
convert a string to a Wiki Word, we split it on punctuation boundaries
and convert each segment to capital case. There are some glitches in this
process. We need to insure the result starts with a letter, and we need
to convert acronyms to words. For example, C<FIGPm> is not valid, but
C<FigPm> is. To get a valid result, we have to convert C<FIG> to C<Fig>.
Also, we have to make sure we don't end a word with C<Tmpl>.

=over 4

=item string

String to convert.

=item suffix

A suffix to add to the string if it has two few humps.

=item RETURN

Incoming string as a WikiWord, with a minimum loss of meaning.

=back

=cut

sub Wikify {
    # Get the parameters.
    my ($string, $suffix) = @_;
    # Declare the return variable.
    my $retVal;
    # Check to see if there's any work to do.
    if (IsWikiWord($string)) {
        $retVal = $string;
    } else {
        # Here we need to convert what we have into a wiki word. First, bust
        # it into segments. Note that underscore (normally a "word" character)
        # is treated as punctuation.
        my @segments = split /[\W_]+/, $string;
        # Insure each segment is capital case.
        for my $segment (@segments) {
            if ($segment =~ /^[A-Z0-9]+$/) {
                # Here we have a probable acronym. Convert it to capital case.
                $segment = ucfirst lc $segment;
            } elsif ($segment =~ /^[a-z0-9]+$/) {
                # This is a lower-case word. Capitalize it.
                $segment = ucfirst $segment;
            }
        }
        # Paste the segments together.
        $retVal = join("", @segments);
        # See if we're done.
        if (! IsWikiWord($retVal) && $suffix) {
            # We're not a WikiWord, but the caller has supplied a suffix.
            # Tack it on.
            $retVal .= $suffix;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 Save

    my $rc = $wiki->Save($title, $web, $category, $content);

Store text for the specified page. If successful, return TRUE. If an
error occurs, return FALSE. In this last case, an error message will be
in the member C<$wiki->{error}>.

=over 4

=item title

Page title, consisting of one or more words. The title will be converted
to the format required by the particular Wiki in use.

=item web

Namespace for this page.

=item category (optional)

Category to contain the page.

=item content

New content for the page.

=item RETURN

Returns TRUE if successful, else FALSE. If FALSE is returned, an error message
will be stashed in the C<error> member.

=back

=cut

sub Save {
    # Get the parameters.
    my ($self, $title, $web, $category, $content) = @_;
    # The first task is to create the real page title. First, we separate the incoming
    # title into words.
    my @words = split /\s+|_+/, $title;
    # Form them into a WikiWord with a web name prefix.
    my $wordTitle = join("", map { ucfirst $_ } @words);
    my $realTitle = "$web.$wordTitle";
    # Set up the parent information (if any).
    my @parentData = ();
    if ($category) {
        push @parentData, parent => $category;
    }
    # Send a request to the wiki server.
    my $ua = $self->{userAgent};
    my $response = $ua->post($self->{url}, { username => $self->{username},
                                             password => $self->{password},
                                             topic    => $realTitle,
                                             text     => $content,
                                             @parentData });
    # Declare the return variable. We assume failure. If we succeed, we'll change
    # the value.
    my $retVal = 0;
    # Save the response content.
    my $message = $response->content;
    Trace("Message returned after Save:\n$message") if T(3);
    # Check for a response error.
    if (! $response->is_success()) {
        # Denote failure, and save the content as the error message.
        if ($message =~ /<p>(.+)<\/p>/si) {
            # Here the message is formatted as HTML. We return the meat.
            $self->{error} = $1;
            $self->{error} =~ s/\n/ /gs;
        } else {
            $self->{error} = $message;
        }
    } elsif ($message =~ /^ERROR:\s*(.*)/) {
        # Here we have an error detected by the plugin script. We get the
        # meat of the message error field.
        $self->{error} = $1;
    } else {
        # Denote success.
        $retVal = 1;
    }
    # Return the result.
    return $retVal;
}


=head3 Bar

    my $line = $wt->Bar;

Return the code for a horizontal bar.

=cut

sub Bar {
    return "---";
}

=head3 HeadParse

    my ($level, $name) = $wiki->HeadParse($line);

Determine whether or not the specified line is a wiki heading. If it is,
return the heading level and the heading name. If it is not, return 0
and C<undef>

=over 4

=item line

Wiki line to parse.

=item RETURN

Returns a two-element list. The first element indicates the heading level, and is
C<0> for a non-heading. The second contains the heading name, or C<undef> for a
non-heading.

=back

=cut

sub HeadParse {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the parameters.
    my ($line) = @_;
    # Assume this is not a heading line.
    my ($level, $name) = (0, undef);
    # Check for the heading format.
    if ($line =~ /^\s*---(\++)\s*(.+)\s*$/) {
        # Here we have a heading line. The heading level is the number of plus
        # signs matched, which is also the length of $1.
        $level = length $1;
        $name = $2;
    }
    # Return the result.
    return ($level, $name);
}

=head3 HeadLevel

    my $level = $wiki->HeadLevel($line);

Return the heading level of the specified line, or 0 if it is not a
heading.

=over 4

=item line

Wiki line to parse.

=item RETURN

Returns C<0> if the line is not a heading, and the heading level otherwise.

=back

=cut

sub HeadLevel {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the parameters.
    my ($line) = @_;
    # Parse the header. We keep the heading level and throw away the text.
    my ($retVal) = HeadParse($line);
    # Return the result.
    return $retVal;
}

=head3 BoldCode

    my $boldCode = $wiki->BoldCode();

Returns the Wiki code for bold text.

=cut

sub BoldCode {
    # Return the result.
    return "*";
}

=head3 ItalicCode

    my $italicCode = $wiki->BoldCode();

Returns the Wiki code for italic text.

=cut

sub ItalicCode {
    # Return the result.
    return "_";
}

=head3 ListCode

    my $listCode = $wiki->ListCode();

Returns the Wiki code for a list element.

=cut

sub ListCode {
    # Return the result.
    return "   * ";
}

=head3 IsWikiWord

    my $flag = WikiTools::IsWikiWord($string);

Return TRUE if the specified string is a [[TWiki.WikiWord][wiki word]], else FALSE.

=over 4

=item string

String to evaluate.

=item RETURN

Returns TRUE if the string conforms to the allowable format for a Wiki page title,
else FALSE.

=back

=cut

sub IsWikiWord {
    # Get the parameters.
    my ($string) = @_;
    # Test the string.
    return $string =~ /^[A-Z]+[a-z]+(?:[A-Z]+[a-zA-Z0-9]*)$/;
}


=head2 Rendering Methods

These are the methods that need to be replicated in any object used for
rendering ERDB documentation.

=head3 Heading

    my $line = $wiki->Heading($level, $text);

Return the code for a heading line at the specified level.

=over 4

=item level

Desired heading level.

=item text

Title for the heading's section.

=item RETURN

Returns a formatted heading line.

=back

=cut

sub Heading {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the parameters.
    my ($level, $text) = @_;
    # Create the heading line.
    my $retVal = "---" . ("+" x $level) . " $text";
    # Return the result.
    return $retVal;
}


=head3 Prolog

    my @lines = $wiki->Prolog();

Returns a set of text lines to put at the beginning of a typical Wiki
output stream.

=cut

sub Prolog {
    # Return the result.
    return ('<noautolink>', '%TOC%');
}

=head3 Epilog

    my @lines = $wiki->Epilog();

Returns a set of text lines to put at the end of a typical Wiki
output stream.

=cut

sub Epilog {
    # Return the result.
    return ('</noautolink>');
}

=head3 Bold

    my $markup = $wiki->Bold($text);

Bold the specified text.

=cut

sub Bold {
    my ($self, $text) = @_;
    return "*$text*";
}

=head3 Italic

    my $markup = $wiki->Italic($text);

Italicize the specified text.

=cut

sub Italic {
    my ($self, $text) = @_;
    return "_" . $text . "_";
}

=head3 LinkMarkup

    my $boldCode = $wiki->LinkMarkup($link, $text);

Returns the Wiki code for a link.

=over 4

=item link

URL or topic name referenced by the link.

=item text (optional)

Text of the link.

=back

=cut

sub LinkMarkup {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the parameters.
    my ($link, $text) = @_;
    # Declare the return variable.
    my $retVal;
    # Check to see if we have text.
    if ($text) {
        # Yes, so we have a two-part link.
        $retVal = "[[$link][$text]]";
    } else {
        # No, so we have a one-part link.
        $retVal = "[[$link]]";
    }
    # Return the result.
    return $retVal;
}

=head3 Table

    my $wikiText = $wiki->Table(@rows);

Create a Wiki table. The parameters are all list references. The first
describes the header row, and the remaining rows are presented
sequentially. This is a very simple table, using only default settings
and with everything left-aligned.

=over 4

=item rows

List of table rows. Each table row is a list reference containing the
cells of the row in column order. The first row is used as the header.

=item RETURN

Returns a string that will generate a Wiki table.

=back

=cut

sub Table {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the parameters.
    my (@rows) = @_;
    # Get the header row.
    my $headers = shift @rows;
    # We put asterisks around the title of each column so that TWiki knows these are headers.
    my $headerRow = "| " . join(" | ", map { "*$_*" } @{$headers}) . " |";
    # Save the header, and build the rest of the rows normally.
    my @rowStrings = $headerRow;
    for my $row (@rows) {
        # Remove line-feeds from the cells.
        my @cells = map { $_ =~ s/\n/ /g; $_ } @{$row};
        # Add them together to make a row.
        push @rowStrings, "| " . join(" | ", @cells) . " |";
    }
    # Put the rows together with blank lines on either side.
    my $retVal = join("\n", "", @rowStrings, "");
    # Return the result.
    return $retVal;
}


=head3 List

    my $wikiText = $wiki->List(@items);

Create a Wiki list. The parameters are all strings that are put into the
list sequentially. The strings are trimmed, and empty entries at the
beginning are deleted. This makes the coding of asides in ERDB a little
more user-friendly.

=over 4

=item items

List of items to be formatted into a wiki list.

=item RETURN

Returns wiki markup text that will display as an unordered list.

=back

=cut

sub List {
    # Allow static calling for backward compatability.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    # Get the list elements, trimmed.
    my (@items) = map { Tracer::Trim($_) } @_;
    # Get the list code.
    my $code = ListCode();
    # Remove any null entries at the beginning.
    while (@items && $items[0] eq "") { shift @items };
    # Format the list.
    my $retVal = join("\n", "", map { "$code$_" } @items);
    # Return the result.
    return $retVal;
}

=head3 Para

    my $markup = $wiki->Para($text);

Create a paragraph from the specified text.

=over 4

=item text

Text to format as a paragraph.

=item RETURN

Returns the text followed by a blank line, so that it is treated as a
paragraph.

=back

=cut

sub Para {
    # Get the parameters.
    my ($self, $text) = @_;
    # Add the blank line.
    my $retVal = "$text\n\n";
    # Return the result.
    return $retVal;
}


=head3 Finalize

    $wiki->Finalize(\@lines);

Finalize a list of lines into a wiki page. This method is not used in the
current implementation, but would be needed by an HTML utility to
generate the table of contents.

=over 4

=item lines

List of lines containing markup. The list is modified in place.

=back

=cut

sub Finalize {
    # Stub.
}


1;
