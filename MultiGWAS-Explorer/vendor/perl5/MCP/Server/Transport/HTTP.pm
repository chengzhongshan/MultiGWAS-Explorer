package MCP::Server::Transport::HTTP;
use Mojo::Base 'MCP::Server::Transport', -signatures;

use Scalar::Util qw(blessed);

sub handle_request ($self, $c) {
  my $request = $c->req->json;
  my $response = $self->server->handle($request, {c => $c, tx => $c->tx});

  return $c->rendered(202) unless defined $response;

  if (blessed($response) && $response->isa('Mojo::Promise')) {
    return $response->then(sub ($resolved) { $c->render(json => $resolved) })
      ->catch(sub ($err) {
        $c->render(
          json => {
            jsonrpc => '2.0',
            error   => {code => -32603, message => "$err"},
          },
          status => 500,
        );
      });
  }

  return $c->render(json => $response);
}

1;

=encoding utf8

=head1 NAME

MCP::Server::Transport::HTTP - HTTP transport for MCP servers

=head1 DESCRIPTION

Small Mojolicious transport used by C<MCP::Server/to_action>. It accepts a JSON-RPC request body and renders the
server response as JSON, returning HTTP 202 for notifications that do not produce a response.

=cut
